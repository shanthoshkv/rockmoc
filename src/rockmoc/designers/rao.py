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
    # Base angles at Lf=0.8: θ_n rises slowly with ε, θ_e falls with ε.
    # A shorter bell turns harder: BOTH θ_n and θ_e grow as Lf drops (Huzel & Huang
    # Fig. 4-16: at ε≈6, 60 % ≈ +6°/+5.5°, 100 % ≈ −4°/−3° relative to 80 %).
    # (v0.9 had θ_e falling with shorter length — inverted; CFD showed the
    # resulting 60 % bell over-expanding and recompressing at the lip.)
    loge = math.log(eps)
    d = Lf - 0.8
    theta_n = 22.0 + 2.5 * loge - 25.0 * d + 25.0 * d * d
    theta_e = 14.0 - 1.8 * loge - 21.25 * d + 31.25 * d * d
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
    length_m: float | None = None,
) -> tuple[Contour, dict]:
    """Design a Rao-type TOP bell contour.

    Parameters
    ----------
    length_fraction:
        Fraction of a 15° conical nozzle length (0.6–1.0). 0.8 is common.
    rc_over_rt:
        Throat wall radius of curvature / Rt.
    length_m:
        Throat-to-exit length overriding ``length_fraction`` × 15° cone length
        (angles still from the clipped ``length_fraction``). Used by TOC.
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
    L = Lf * L_cone if length_m is None else float(length_m)

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

    # Rao TOP parabola N→E: quadratic Bézier whose control point Q is the
    # intersection of the tangents at N (θ_n) and E (θ_e). End slopes are exact,
    # with no overshoot (the v0.9 cubic Hermite overshot θ_e by ~1°).
    n_par = max(30, n_points - n_arc)
    t = np.linspace(0.0, 1.0, n_par).reshape(-1, 1)
    tn, te = math.tan(theta_n), math.tan(theta_e)
    xq = (r_e - r_n + x_n * tn - x_e * te) / (tn - te)
    P0, P2 = np.array([x_n, r_n]), np.array([x_e, r_e])
    shape = "quadratic_bezier"
    if x_n < xq < x_e:
        Q = np.array([xq, r_n + (xq - x_n) * tn])
        P = (1 - t) ** 2 * P0 + 2 * (1 - t) * t * Q + t**2 * P2
        dP = 2 * (1 - t) * (Q - P0) + 2 * t * (P2 - Q)
    else:
        # ponytail: very short / low-ε bells admit no TOP parabola for these θ_n, θ_e
        # (tangents cross outside N–E); cubic Hermite through the same end slopes instead.
        shape = "cubic_hermite_fallback"
        L_ne = x_e - x_n
        m0 = (
            0.4
            * L_ne
            / max(math.cos(theta_n), 0.2)
            * np.array([math.cos(theta_n), math.sin(theta_n)])
        )
        m1 = (
            0.4
            * L_ne
            / max(math.cos(theta_e), 0.2)
            * np.array([math.cos(theta_e), math.sin(theta_e)])
        )
        P = (
            (2 * t**3 - 3 * t**2 + 1) * P0
            + (t**3 - 2 * t**2 + t) * m0
            + (-2 * t**3 + 3 * t**2) * P2
            + (t**3 - t**2) * m1
        )
        dP = (
            (6 * t**2 - 6 * t) * P0
            + (3 * t**2 - 4 * t + 1) * m0
            + (-6 * t**2 + 6 * t) * P2
            + (3 * t**2 - 2 * t) * m1
        )
    x_par, r_par = P[:, 0], P[:, 1]

    # Concatenate (skip duplicate junction)
    x = np.concatenate([x_arc[:-1], x_par])
    r = np.concatenate([r_arc[:-1], r_par])
    x[-1], r[-1] = x_e, r_e
    r = np.maximum.accumulate(r)  # non-decreasing radius (Hermite fallback can dip)
    # analytic wall angle (arc: φ, parabola: Bézier tangent)
    th = np.concatenate([phi[:-1], np.arctan2(dP[:, 1], dP[:, 0])])
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
        "parabola": shape,
        "reference": "Rao TOP engineering correlations (Huzel & Huang trends) / quadratic Bezier parabola",
    }
    return contour, meta


def design_rao(
    spec: DesignSpec, *, length_fraction: float = 0.80, store_net: bool = False
) -> DesignResult:
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
