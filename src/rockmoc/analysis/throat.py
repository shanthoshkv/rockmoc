"""Transonic throat initial-value line (IVL) from Sauer's solution.

The velocity field near a circular-arc throat is Sauer's small-perturbation
solution (see :mod:`rockmoc.moc.transonic` for equations and references).
The IVL is the vertical line where the axis Mach equals ``axis_mach``; it is
weakly supersonic everywhere and spacelike, so MOC can march from it.
Sauer is accurate for Rc/Rt ≳ 2. This is **not** CFD.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from rockmoc.exceptions import DesignError
from rockmoc.moc.transonic import sauer_alpha, sauer_line, sauer_offset
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


def sauer_throat_ivl(
    gas: PerfectGas,
    *,
    throat_radius: float,
    rc_over_rt: float = 1.5,
    n_points: int = 41,
    axis_mach: float = 1.02,
    delta: float = 1.0,
) -> ThroatIVL:
    """Sauer IVL on the vertical line where the axis Mach is ``axis_mach``."""
    if throat_radius <= 0 or rc_over_rt <= 0:
        raise DesignError("throat_radius and rc_over_rt must be positive")
    if n_points < 5:
        raise DesignError("n_points must be ≥ 5")
    if axis_mach <= 1.0:
        raise DesignError("axis_mach must be > 1")

    Rt = float(throat_radius)
    Rc = rc_over_rt * Rt
    g = float(gas.gamma)
    try:
        x0, phi_w, line = sauer_line(g, Rt, Rc, axis_mach=axis_mach, n=n_points - 1, delta=delta)
    except ValueError as exc:
        raise DesignError(str(exc)) from exc
    r = np.array([p[0] for p in line])
    M = np.array([p[1] for p in line])
    theta = np.array([p[2] for p in line])
    points = [
        FlowPoint.from_mach_theta(
            x0, float(r[i]), float(M[i]), float(theta[i]), gas,
            kind=PointKind.THROAT if i == n_points - 1 else PointKind.INTERIOR,
        )
        for i in range(n_points)
    ]
    return ThroatIVL(
        x=np.full(n_points, x0),
        r=r,
        M=M,
        theta=theta,
        points=points,
        meta={
            "model": "sauer",
            "Rc_over_Rt": rc_over_rt,
            "x_over_Rt": x0 / Rt,
            "sonic_point_offset_over_Rt": sauer_offset(g, Rt, Rc, delta) / Rt,
            "alpha_per_m": sauer_alpha(g, Rt, Rc, delta),
            "wall_angle_deg": math.degrees(phi_w),
            "gamma": g,
            "references": "Sauer 1947; Mendelson NASA TM X-53084 Eq. 2.6",
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
