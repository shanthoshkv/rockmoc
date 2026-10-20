"""Rao thrust-optimized parabolic (TOP) / bell approximation.

Reference: G.V.R. Rao, exhaust nozzle contour papers (1958–1961).
Industry practice: circular throat arc + parabola to exit with design
wall angles θ_n (inflection) and θ_e (exit) correlated to area ratio and
percent length of a 15° conical nozzle.

This is an engineering TOP model using widely used Huzel & Huang / Rao
angle correlations — not a full method-of-characteristics optimization.
"""

from __future__ import annotations

import math

import numpy as np

from rockmoc.designers.base import DesignResult, DesignSpec
from rockmoc.exceptions import DesignError
from rockmoc.geometry.contour import Contour
from rockmoc.performance.metrics import compute_metrics
from rockmoc.thermo.perfect_gas import PerfectGas


def _rao_angles_deg(area_ratio: float, length_fraction: float) -> tuple[float, float]:
    """Approximate θ_n, θ_e [deg] vs ε and % length of 15° cone.

    Correlations fitted to published Rao TOP charts (engineering fit):
    length_fraction in (0.6, 1.0) typical; 0.8 common.
    """
    eps = max(float(area_ratio), 1.5)
    Lf = float(np.clip(length_fraction, 0.60, 1.0))
    # Base angles at Lf=0.8
    # θ_n decreases slowly with ε; θ_e decreases with ε
    loge = math.log(eps)
    theta_n = 22.0 + 2.5 * loge - 12.0 * (Lf - 0.8)  # ~20–30 deg
    theta_e = 14.0 - 1.8 * loge + 10.0 * (Lf - 0.8)  # ~5–15 deg
    theta_n = float(np.clip(theta_n, 15.0, 50.0))
    theta_e = float(np.clip(theta_e, 2.0, 25.0))
    if theta_e > theta_n - 3.0:
        theta_e = theta_n - 3.0
    return theta_n, theta_e


def conical_length(throat_radius: float, area_ratio: float, cone_angle_deg: float = 15.0) -> float:
    """Axial length of a conical nozzle from throat to exit."""
    Re = throat_radius * math.sqrt(area_ratio)
    return (Re - throat_radius) / math.tan(math.radians(cone_angle_deg))


def design_rao_top(
    gas: PerfectGas,
    *,
    throat_radius: float,
    exit_mach: float | None = None,
    area_ratio: float | None = None,
    length_fraction: float = 0.80,
    rc_over_rt: float = 1.5,
    n_points: int = 200,
) -> tuple[Contour, dict]:
    """Design a Rao-type TOP bell contour.

    Parameters
    ----------
    length_fraction:
        Fraction of a 15° conical nozzle length (0.6–1.0). 0.8 is common.
    rc_over_rt:
        Throat wall radius of curvature / Rt.
    """
    if throat_radius <= 0:
        raise DesignError("throat_radius must be positive")
    if exit_mach is None and area_ratio is None:
        raise DesignError("provide exit_mach or area_ratio")
    if exit_mach is not None:
        eps = float(gas.area_ratio(exit_mach))
        Me = float(exit_mach)
    else:
        eps = float(area_ratio)
        Me = float(gas.mach_from_area_ratio(eps, branch="supersonic"))

    Rt = float(throat_radius)
    Re = Rt * math.sqrt(eps)
    Rc = rc_over_rt * Rt
    Lf = float(np.clip(length_fraction, 0.60, 1.0))
    L_cone = conical_length(Rt, eps, 15.0)
    L = Lf * L_cone

    theta_n_deg, theta_e_deg = _rao_angles_deg(eps, Lf)
    theta_n = math.radians(theta_n_deg)
    theta_e = math.radians(theta_e_deg)

    # Entrance: circular arc from throat (θ=0) to θ_n
    # Arc: x = Rc sinφ, r = Rt + Rc(1-cosφ), φ: 0 → θ_n
    n_arc = max(20, n_points // 5)
    phi = np.linspace(0.0, theta_n, n_arc)
    x_arc = Rc * np.sin(phi)
    r_arc = Rt + Rc * (1.0 - np.cos(phi))

    # Junction N (end of arc)
    x_n = float(x_arc[-1])
    r_n = float(r_arc[-1])

    # Exit E
    x_e = L
    r_e = Re
    if x_e <= x_n:
        x_e = x_n + 0.5 * Rt

    # Parabola from N to E with end slopes tan(θ_n), tan(θ_e)
    # Use cubic Hermite (Bezier-like) which matches Rao parabola family well:
    # r(t), x(t) with t in [0,1]
    n_par = max(30, n_points - n_arc)
    t = np.linspace(0.0, 1.0, n_par)
    # Control length along tangents
    c1 = 0.40 * (x_e - x_n) / max(math.cos(theta_n), 0.2)
    c2 = 0.40 * (x_e - x_n) / max(math.cos(theta_e), 0.2)
    # Cubic Hermite:
    # P(t) = (2t^3-3t^2+1)P0 + (t^3-2t^2+t)m0 + (-2t^3+3t^2)P1 + (t^3-t^2)m1
    P0 = np.array([x_n, r_n])
    P1 = np.array([x_e, r_e])
    m0 = c1 * np.array([math.cos(theta_n), math.sin(theta_n)])
    m1 = c2 * np.array([math.cos(theta_e), math.sin(theta_e)])
    t = t.reshape(-1, 1)
    h00 = 2 * t**3 - 3 * t**2 + 1
    h10 = t**3 - 2 * t**2 + t
    h01 = -2 * t**3 + 3 * t**2
    h11 = t**3 - t**2
    P = h00 * P0 + h10 * m0 + h01 * P1 + h11 * m1
    x_par = P[:, 0]
    r_par = P[:, 1]

    # Concatenate (skip duplicate junction)
    x = np.concatenate([x_arc[:-1], x_par])
    r = np.concatenate([r_arc[:-1], r_par])
    # Ensure end exact
    x[-1], r[-1] = x_e, r_e
    r = np.maximum.accumulate(r)  # non-decreasing radius
    # Wall angle
    th = np.arctan2(np.gradient(r), np.gradient(x))
    # Mach estimate along wall from area
    M = np.ones_like(x)
    for i in range(len(x)):
        ar = max((r[i] / Rt) ** 2, 1.0)
        try:
            M[i] = float(gas.mach_from_area_ratio(ar, branch="supersonic"))
        except Exception:
            M[i] = 1.0 + (Me - 1.0) * i / max(len(x) - 1, 1)

    contour = Contour(x=x, r=r, theta=th, M=M)
    meta = {
        "method": "rao_top_parabola",
        "length_fraction": Lf,
        "L_cone_15deg_m": L_cone,
        "length_m": L,
        "theta_n_deg": theta_n_deg,
        "theta_e_deg": theta_e_deg,
        "area_ratio": eps,
        "exit_mach": Me,
        "rc_over_rt": rc_over_rt,
        "reference": "Rao TOP engineering correlations / cubic Hermite bell",
    }
    return contour, meta


def design_rao(spec: DesignSpec, *, length_fraction: float = 0.80, store_net: bool = False) -> DesignResult:
    """DesignSpec-compatible Rao TOP designer."""
    spec.validate()
    Me = spec.resolved_exit_mach()
    contour, meta = design_rao_top(
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach=Me,
        length_fraction=length_fraction,
        rc_over_rt=spec.throat_rc_over_rt,
        n_points=max(100, spec.n_characteristics * 2),
    )
    metrics = compute_metrics(
        contour,
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach_design=Me,
    )
    return DesignResult(
        kind="rao",
        spec=spec,
        exit_mach_design=Me,
        contour=contour,
        metrics=metrics,
        net=None,
        meta={**meta, "designer": "rao"},
    )
