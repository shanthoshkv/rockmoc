"""Transonic throat IVL — Sauer / Hall / Kliegel first-order theory.

References (see ``refs/``)
--------------------------
* Sauer — small-disturbance throat solution
* Hall — higher-order axisymmetric series
* Kliegel & Levine (1969) — confirms first-order throat conditions match
  Sauer/Hall; useful for small throat radius of curvature
* Dutton (1980) DTIC — review of Hall-type series for nozzle throats

First-order axisymmetric form (engineering)
-------------------------------------------
With ε = R_t / R_c (throat radius / wall curvature radius) and local sonic
speed as reference, the leading-order expansion of the velocity field near
the geometric throat gives a weakly supersonic data line for MOC start::

    M²(x, r) − 1 ≈ α · (x / R_t) + β · (r / R_t)²

where α, β depend on γ and ε (Hall factor improves α for moderate ε).

This is **not** CFD — it seeds MOC downstream of the sonic line.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from rockmoc.exceptions import DesignError
from rockmoc.moc.types import FlowPoint, PointKind
from rockmoc.thermo.perfect_gas import PerfectGas


@dataclass(slots=True)
class ThroatIVL:
    """Initial-value line just past the geometric throat."""

    x: NDArray[np.floating]
    r: NDArray[np.floating]
    M: NDArray[np.floating]
    theta: NDArray[np.floating]
    points: list[FlowPoint]
    meta: dict

    def summary(self) -> dict:
        return {
            "n_points": int(len(self.x)),
            "M_axis": float(self.M[0]),
            "M_wall": float(self.M[-1]),
            "x_mean": float(np.mean(self.x)),
            **self.meta,
        }


def _sauer_hall_coefficients(gamma: float, eps: float) -> tuple[float, float, float]:
    """Return (α, β, hall) for M²−1 = α(x/Rt) + β(r/Rt)².

    First-order Sauer/Kliegel axisymmetric grouping with Hall correction
    on the axial coefficient for finite ε = Rt/Rc.
    """
    g = float(gamma)
    # Leading axial acceleration (order √ε)
    # Classical: α ~ √((γ+1) ε)  (dimensionless x/Rt scaling)
    hall = 1.0 + 0.5 * eps + 0.15 * eps**2  # mild Hall-style series
    alpha = math.sqrt(max((g + 1.0) * eps, 1e-16)) * hall
    # Radial: from continuity / irrotational constraint at first order
    # β = α * (γ+1)/4   is a common Hall/Sauer radial coupling sign/magnitude
    # Positive β ⇒ M increases slightly toward wall on a fixed-x plane
    beta = alpha * (g + 1.0) / 4.0
    return alpha, beta, hall


def sauer_throat_ivl(
    gas: PerfectGas,
    *,
    throat_radius: float,
    rc_over_rt: float = 1.5,
    n_points: int = 41,
    x_over_rt: float = 0.05,
    min_supersonic: float = 0.02,
) -> ThroatIVL:
    """Build a Sauer/Hall/Kliegel slightly-supersonic IVL on x = const.

    Parameters
    ----------
    min_supersonic:
        Minimum (M²−1) enforced on the IVL so MOC stays hyperbolic.
    """
    if throat_radius <= 0 or rc_over_rt <= 0:
        raise DesignError("throat_radius and rc_over_rt must be positive")
    if n_points < 5:
        raise DesignError("n_points must be ≥ 5")

    Rt = float(throat_radius)
    Rc = rc_over_rt * Rt
    g = gas.gamma
    eps = Rt / Rc
    alpha, beta, hall = _sauer_hall_coefficients(g, eps)

    x0 = x_over_rt * Rt
    # Ensure weakly supersonic everywhere on the plane
    for _ in range(12):
        m2m1_axis = alpha * (x0 / Rt)
        m2m1_wall = alpha * (x0 / Rt) + beta * 1.0
        if min(m2m1_axis, m2m1_wall) >= min_supersonic:
            break
        x0 *= 1.35
    x0 = max(x0, 0.02 * Rt)

    r_wall = Rt  # first-order: cylindrical section at geometric throat radius
    r = np.linspace(0.0, r_wall, n_points)
    x = np.full(n_points, x0)
    M2m1 = alpha * (x0 / Rt) + beta * (r / Rt) ** 2
    M2m1 = np.maximum(M2m1, min_supersonic)
    M = np.sqrt(1.0 + M2m1)

    # Flow angle: from irrotationality / wall slope consistency
    # At wall, circular arc past throat: φ ≈ arcsin(x/Rc) ≈ x/Rc for small x
    phi_w = math.asin(min(1.0, x0 / max(Rc, 1e-12)))
    # Linear θ(r) from axis (0) to wall (φ_w) — leading-order Sauer
    theta = phi_w * (r / max(r_wall, 1e-12))

    points: list[FlowPoint] = []
    for i in range(n_points):
        points.append(
            FlowPoint.from_mach_theta(
                float(x[i]),
                float(r[i]),
                float(M[i]),
                float(theta[i]),
                gas,
                kind=PointKind.THROAT if i == n_points - 1 else PointKind.INTERIOR,
            )
        )

    return ThroatIVL(
        x=x,
        r=r,
        M=M,
        theta=theta,
        points=points,
        meta={
            "model": "sauer_hall_kliegel_first_order",
            "Rc_over_Rt": rc_over_rt,
            "eps_Rt_over_Rc": eps,
            "x_over_Rt": float(x0 / Rt),
            "gamma": g,
            "alpha": alpha,
            "beta": beta,
            "hall_factor": hall,
            "references": "Sauer; Hall; Kliegel & Levine 1969",
        },
    )


def circular_arc_throat_wall(
    gas: PerfectGas,
    *,
    throat_radius: float,
    rc_over_rt: float = 1.5,
    theta_end_deg: float = 30.0,
    n_points: int = 40,
    throat_mach: float = 1.05,
) -> list[FlowPoint]:
    """Circular-arc divergent wall from geometric throat to wall angle θ_end.

    Geometry (standard rocket nozzle entrance)::

        x = Rc · sin φ
        r = Rt + Rc · (1 − cos φ)
        wall angle = φ
    """
    Rt = throat_radius
    Rc = rc_over_rt * Rt
    theta_end = math.radians(theta_end_deg)
    pts: list[FlowPoint] = []
    for i in range(n_points + 1):
        phi = theta_end * i / n_points
        x = Rc * math.sin(phi)
        r = Rt + Rc * (1.0 - math.cos(phi))
        ar = max((r / Rt) ** 2, 1.0)
        try:
            M = float(gas.mach_from_area_ratio(ar, branch="supersonic"))
        except Exception:
            M = throat_mach
        M = max(M, throat_mach if i == 0 else M)
        if i == 0:
            M = max(throat_mach, 1.0 + 1e-8)
        pts.append(
            FlowPoint.from_mach_theta(
                x, r, M, phi, gas, kind=PointKind.THROAT if i == 0 else PointKind.WALL
            )
        )
    return pts
