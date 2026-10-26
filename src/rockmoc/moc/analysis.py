"""Direct MOC analysis of a *given* wall, with internal-shock detection.

The designers build a wall from a characteristic net; this module does the
inverse: march the same smooth-throat lattice (transonic start line, interior,
axis and wall unit processes) against an arbitrary wall polyline, e.g. a TOP
parabola or a cone. Wall-reflected compression waves that are stronger than
the incoming expansion make neighbouring characteristics of one family cross;
the first crossing (a sign change of the cell Jacobian ∂(x, r)/∂(α, β)) marks
where an internal shock forms (Frey et al., EUCASS 2017: "converging
characteristics indicate an internal shock"). This is a detector, not shock
fitting: downstream of the first crossing the isentropic solution is invalid.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import make_interp_spline

from rockmoc.exceptions import DesignError, MeshError, ThermoError
from rockmoc.geometry.contour import Contour
from rockmoc.moc._common import _UNIT_ITER, _UNIT_TOL, _fp
from rockmoc.moc.equations import source_plus
from rockmoc.moc.mesh import axis_point, interior_point
from rockmoc.moc.transonic import HALL_MIN_RC_RT, KL_MIN_RC_RT, hall_line, kl_line, sauer_line
from rockmoc.moc.types import FlowPoint, PointKind
from rockmoc.thermo.gas_protocol import MOCGas


@dataclass
class WallAnalysis:
    has_internal_shock: bool
    shock_origin: tuple[float, float] | None  # (x, r) of the first characteristic crossing [m]
    wall: list[FlowPoint]
    points: dict = field(repr=False)
    rc_over_rt: float = math.nan
    meta: dict = field(default_factory=dict)


def _throat_rc(x: np.ndarray, r: np.ndarray, i_t: int, Rt: float) -> float:
    """Downstream throat radius of curvature from r − Rt ≈ x²/(2 Rc)."""
    xs, rs = x[i_t + 1:] - x[i_t], r[i_t + 1:] - r[i_t]
    m = (xs > 0) & (xs < 0.25 * Rt) & (rs > 0)
    if m.sum() < 2:
        raise DesignError("cannot estimate the throat radius of curvature (refine the contour near the throat)")
    return float(np.median(xs[m] ** 2 / (2.0 * rs[m])))


def analyze_wall(
    contour: Contour,
    gas: MOCGas,
    *,
    n_chars: int = 40,
    throat_mach: float = 1.02,
    rc_over_rt: float | None = None,
    delta: float = 1.0,
    axis_band: float = 0.1,
) -> WallAnalysis:
    """March characteristics against ``contour`` and flag internal shocks.

    ``rc_over_rt`` defaults to a fit of the given points just downstream of the
    throat; pass it when the contour is coarse there. Crossings within
    ``axis_band``·Rt of the axis are ignored (numerical, see module notes).
    """
    xc, rc = np.asarray(contour.x, dtype=np.float64), np.asarray(contour.r, dtype=np.float64)
    i_c = int(np.argmin(rc))
    # Rc from the given points (resampled points would linearise a coarse arc)
    Rc = _throat_rc(xc, rc, i_c, float(rc[i_c])) if rc_over_rt is None else rc_over_rt * float(rc[i_c])
    # C² cubic spline through the given points: a linear polyline's slope kinks would
    # each launch a spurious wave and corrupt the crossing test
    keep = np.concatenate([[True], np.diff(xc) > 1e-12 * (1.0 + np.abs(xc[1:]))])
    spl = make_interp_spline(xc[keep], rc[keep], k=3)
    x0 = np.linspace(xc[0], xc[-1], max(4000, 4 * len(xc)))
    r0 = spl(x0)
    i_t = int(np.argmin(r0))
    Rt = float(r0[i_t])
    xw, rw = x0[i_t:] - x0[i_t], r0[i_t:]
    # Wall angle from the contour's own θ (exact for designer walls) rather than
    # from differentiating r(x), which is biased near a coarsely sampled arc.
    th_c = np.asarray(contour.theta, dtype=np.float64)[keep]
    thw = make_interp_spline(xc[keep], th_c, k=1)(x0[i_t:])
    m = max(int(n_chars), 6)

    R = Rc / Rt
    fns = [sauer_line]
    if abs(delta) >= 0.5:
        fns = ([hall_line] if R >= HALL_MIN_RC_RT else [kl_line] if R >= KL_MIN_RC_RT else []) + fns
    g_t = float(gas.gamma_at_mach(1.0))
    M0 = max(throat_mach, 1.0 + 1e-6)
    line_fn, line = None, None
    for fn in fns:
        try:
            x_line, phi, cand = fn(g_t, Rt, Rc, axis_mach=M0, n=m, delta=delta)
        except ValueError:
            continue
        # the series start line assumes the circular arc: reject it past the arc tangency
        if abs(phi - float(np.interp(x_line, xw, thw))) < math.radians(2.0):
            line_fn, line = fn, cand
            break
    if line is None:  # uniform near-sonic line on the throat plane (as the designers' fallback)
        x_line, line_fn = 0.0, None
        line = [(Rt * j / m, M0, 0.0) for j in range(m + 1)]
    P: dict[tuple[int, int], FlowPoint] = {}
    for j, (rj, Mj, th) in enumerate(line):
        kind = PointKind.AXIS if j == 0 else (PointKind.WALL if j == m else PointKind.INTERIOR)
        P[(m - j, j)] = _fp(x_line, rj, th, float(gas.prandtl_meyer(Mj)), gas, kind)
    x_end = float(xw[-1])

    def wall_hit(b: FlowPoint, al: float) -> tuple[float, float]:
        """First intersection of the ray from b at angle al with the wall polyline."""
        dx, dr = math.cos(al), math.sin(al)
        k0 = max(int(np.searchsorted(xw, b.x)) - 1, 0)
        for k in range(k0, len(xw) - 1):
            ex, er = xw[k + 1] - xw[k], rw[k + 1] - rw[k]
            den = dx * er - dr * ex
            if abs(den) < 1e-300:
                continue
            t = ((xw[k] - b.x) * er - (rw[k] - b.r) * ex) / den
            u = ((xw[k] - b.x) * dr - (rw[k] - b.r) * dx) / den
            if t > 0.0 and -1e-12 <= u <= 1.0 + 1e-12:
                return b.x + t * dx, b.r + t * dr
        raise MeshError("C+ leaves the nozzle through the exit")

    def wall_point(b: FlowPoint) -> FlowPoint:
        al = b.theta + b.mu
        w = None
        for _ in range(_UNIT_ITER):
            x, r = wall_hit(b, al)
            th = float(np.interp(x, xw, thw))
            mu_p = b.mu if w is None else 0.5 * (b.mu + w.mu)
            kp = (b.theta - b.nu) + source_plus(0.5 * (b.theta + th), mu_p, 0.5 * (b.r + r), delta=delta) * (x - b.x)
            w_new = _fp(x, r, th, max(th - kp, 0.0), gas, PointKind.WALL)
            if w is not None and abs(w_new.nu - w.nu) < _UNIT_TOL:
                return w_new
            w = w_new
            al = 0.5 * (b.theta + b.mu + w.theta + w.mu)
        return w  # type: ignore[return-value]

    def point(a: int, b: int) -> FlowPoint:
        if b - a == m:
            return wall_point(P[(a, b - 1)])
        if a - b == m:
            return axis_point(P[(a - 1, b)], gas, delta=delta)
        return interior_point(P[(a - 1, b)], P[(a, b - 1)], gas, delta=delta)

    failures = 0
    for s_ in range(m + 1, m + 1 + 400 * m):
        made = 0
        for a in range(max((s_ - m + 1) // 2, 0), (s_ + m) // 2 + 1):
            b = s_ - a
            if abs(a - b) > m or (a, b) in P:
                continue
            try:
                p = point(a, b)
            except KeyError:  # upstream point outside the nozzle
                continue
            except (MeshError, ThermoError, ValueError):
                # ponytail: near-axis sinθ/r overflow → drop the source term for this one
                # interior point (planar update) rather than stall the whole march
                try:
                    if abs(a - b) == m:
                        raise MeshError("boundary point")
                    up = P[(a - 1, b)]
                    p = interior_point(up, P[(a, b - 1)], gas, delta=0.0)
                    failures += 1
                    if not (math.isfinite(p.x) and math.isfinite(p.r) and p.x >= up.x and p.r >= 0.0):
                        continue
                except (MeshError, ThermoError, ValueError, KeyError):
                    failures += 1
                    continue
            if p.x <= x_end * (1.0 + 1e-9):
                P[(a, b)] = p
                made += 1
        if made == 0:
            break

    # Cell orientation: cross((P(a+1,b) − P(a,b)), (P(a,b+1) − P(a,b))) keeps one sign in a
    # shock-free net; a flip means two characteristics of one family have crossed.
    signs, cells = [], []
    for (a, b), p in P.items():
        q1, q2 = P.get((a + 1, b)), P.get((a, b + 1))
        if q1 is None or q2 is None:
            continue
        cr = (q1.x - p.x) * (q2.r - p.r) - (q1.r - p.r) * (q2.x - p.x)
        signs.append(cr)
        cells.append(p)
    signs_a = np.array(signs)
    ref = np.sign(np.median(signs_a)) if len(signs_a) else 1.0
    # near-axis cells are excluded: the sinθ/r terms there amplify wall-resolution noise
    bad = [c for c, sg in zip(cells, signs_a) if sg * ref <= 0.0 and c.r > axis_band * Rt]
    # ponytail: a shock is a coherent cluster of crossings; isolated flips (≲0.2 % of cells on
    # shock-free ideal walls) are wall-resolution noise. 1 % threshold; shock fitting would replace it.
    origin = min(((c.x, c.r) for c in bad), default=None) if len(bad) >= 0.01 * max(len(cells), 1) else None
    wall = sorted((p for (a, b), p in P.items() if b - a == m), key=lambda p: p.x)
    return WallAnalysis(
        has_internal_shock=origin is not None,
        shock_origin=(origin[0] + float(x0[i_t]), origin[1]) if origin else None,
        wall=wall,
        points=P,
        rc_over_rt=R,
        meta={"n_points": len(P), "unit_failures": failures, "n_crossed_cells": len(bad),
              "initial_line": line_fn.__name__ if line_fn else "uniform", "wall_coverage": (wall[-1].x / x_end) if wall else 0.0},
    )
