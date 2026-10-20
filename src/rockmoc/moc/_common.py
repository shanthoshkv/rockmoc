"""Shared MOC helpers: point construction, line intersection, mass flow, contour utilities."""

from __future__ import annotations

import functools
import math

import numpy as np

from rockmoc.exceptions import collect_fallbacks, warn_fallback
from rockmoc.geometry.contour import Contour
from rockmoc.moc.types import CharNet, FlowPoint, PointKind
from rockmoc.thermo.gas_protocol import MOCGas


def _pm_state(nu: float, gas: MOCGas) -> tuple[float, float, float]:
    """(M, μ, ν_clamped) from Prandtl–Meyer angle."""
    nu_c = max(0.0, float(nu))
    nu_max = float(gas.prandtl_meyer_max()) * 0.999
    if nu_c > nu_max:
        warn_fallback("Prandtl–Meyer angle clamped to ν_max (over-expansion in the net)")
        nu_c = nu_max
    if nu_c <= 1e-16:
        M = 1.01
        nu_c = float(gas.prandtl_meyer(M))
    else:
        try:
            M = float(gas.mach_from_prandtl_meyer(nu_c))
        except Exception:
            warn_fallback("inverse Prandtl–Meyer failed near ν_max; approximate Mach used")
            M = 1.01 + 20.0 * (nu_c / max(nu_max, 1e-6))
            M = min(M, 50.0)
            nu_c = float(gas.prandtl_meyer(M))
    mu = math.asin(min(1.0, 1.0 / max(M, 1.0 + 1e-12)))
    return M, mu, nu_c


def _fp(
    x: float,
    r: float,
    theta: float,
    nu: float,
    gas: MOCGas,
    kind: PointKind,
) -> FlowPoint:
    M, mu, nu_c = _pm_state(nu, gas)
    g_loc = float(gas.gamma_at_mach(M))
    return FlowPoint(
        x=float(x),
        r=max(float(r), 0.0),
        theta=float(theta),
        nu=nu_c,
        M=M,
        mu=mu,
        kind=kind,
        gamma_local=g_loc,
    )


def _find_xy(
    xy_top: tuple[float, float],
    xy_bot: tuple[float, float],
    c_neg: float,
    c_pos: float,
) -> tuple[float, float]:
    """Intersection of C− through top (angle c_neg) and C+ through bot (c_pos).

    Angles are characteristic inclination angles [rad], not slopes.
    Formula matches standard planar MOC references (nategphillips/moc geometry).
    """
    ta = math.tan(c_neg)
    tb = math.tan(c_pos)
    den = ta - tb
    if abs(den) < 1e-14:
        return xy_bot[0] + 1e-9, max(xy_bot[1], 0.0)
    x = (xy_top[0] * ta - xy_bot[0] * tb + xy_bot[1] - xy_top[1]) / den
    y = (ta * tb * (xy_top[0] - xy_bot[0]) + ta * xy_bot[1] - tb * xy_top[1]) / den
    return float(x), max(float(y), 0.0)


_UNIT_ITER = 25


_UNIT_TOL = 1e-11


def _poly_mass(pts: list[FlowPoint], gas: MOCGas, *, planar: bool = False) -> list[float]:
    """Cumulative mass flow ∫ρV·n dA / (ρ*a*) along a polyline (per unit depth if planar)."""
    out = [0.0]
    for p, q in zip(pts[:-1], pts[1:]):
        th = 0.5 * (p.theta + q.theta)
        ar = float(gas.area_ratio(0.5 * (p.M + q.M)))
        dn = math.cos(th) * (q.r - p.r) - math.sin(th) * (q.x - p.x)
        w = 1.0 if planar else math.pi * (p.r + q.r)
        out.append(out[-1] + w * dn / ar)
    return out


def _wall_to_contour(wall: list[FlowPoint], Rt: float, *, axial_exit: bool = True) -> Contour:
    xs = np.array([p.x for p in wall], dtype=np.float64)
    rs = np.array([p.r for p in wall], dtype=np.float64)
    th = np.array([p.theta for p in wall], dtype=np.float64)
    Mv = np.array([p.M for p in wall], dtype=np.float64)
    rs[0] = Rt
    xs[0] = 0.0
    th[0] = max(th[0], th[1] if len(th) > 1 else th[0])  # throat after turn
    for i in range(1, len(xs)):
        if xs[i] <= xs[i - 1]:
            xs[i] = xs[i - 1] + 1e-14
        if rs[i] < rs[i - 1] * 0.999:
            rs[i] = rs[i - 1]
    if axial_exit:
        th[-1] = 0.0
    return Contour(x=xs, r=rs, theta=th, M=Mv)


def _isotropic_match_exit(contour: Contour, Rt: float, r_exit_target: float) -> Contour:
    if contour.r[-1] <= Rt + 1e-15 or r_exit_target <= Rt:
        return contour
    scale = (r_exit_target - Rt) / (float(contour.r[-1]) - Rt)
    x0 = float(contour.x[0])
    x = x0 + scale * (contour.x - x0)
    r = Rt + scale * (contour.r - Rt)
    r[-1] = r_exit_target
    return Contour(x=x, r=r, theta=contour.theta, M=contour.M)


def _scale_net(net: CharNet, Rt: float, scale: float, x0: float = 0.0) -> None:
    if abs(scale - 1.0) < 1e-14:
        return

    def _sc(p: FlowPoint) -> FlowPoint:
        return p.copy(
            x=x0 + scale * (p.x - x0),
            r=(Rt + scale * (p.r - Rt)) if p.r > 1e-14 else 0.0,
        )

    net.points = [_sc(p) for p in net.points]
    net.axis = [_sc(p) for p in net.axis]
    net.wall = [_sc(p) for p in net.wall]
    net.fan = [_sc(p) for p in net.fan]
    net.interior = [[_sc(p) for p in row] for row in net.interior]


def blend_contours(c_lo: Contour, c_hi: Contour, w: float, *, n: int = 200, axial_exit: bool = True) -> Contour:
    """Arc-length-parametric blend (1 − w)·c_lo + w·c_hi of two nearby designs.

    Used where a design parameter is discrete on the characteristic lattice
    (arc attachment node): the two bracketing designs are interpolated.
    """
    a, b = c_lo.resample(n), c_hi.resample(n)
    th = (1.0 - w) * a.theta + w * b.theta
    th[0] = 0.0
    if axial_exit:
        th[-1] = 0.0
    return Contour(x=(1.0 - w) * a.x + w * b.x, r=(1.0 - w) * a.r + w * b.r, theta=th,
                   M=(1.0 - w) * a.M + w * b.M)


def _records_fallbacks(fn):  # type: ignore[no-untyped-def]
    """Run a designer, collecting MOCFallbackWarning messages into result meta['warnings']."""

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):  # type: ignore[no-untyped-def]
        with collect_fallbacks() as msgs:
            res = fn(*args, **kwargs)
        meta = res.meta if hasattr(res, "meta") else res[-1]
        meta.setdefault("warnings", [])
        meta["warnings"] += [m for m in msgs if m not in meta["warnings"]]
        return res

    return wrapper
