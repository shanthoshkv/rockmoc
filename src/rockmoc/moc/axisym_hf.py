"""High-fidelity MOC nozzle design — full characteristic net (kernel + wall).

Algorithm sources (methods, not code copies)
--------------------------------------------
* Anderson, *Modern Compressible Flow* — planar MLN construction
* Zucrow & Hoffman — axisymmetric geometric source terms
* NASA TM-X-1502 — planar sharp-throat rules (θ_max = ν_e/2)
* Planar mesh topology consistent with open references:
  nategphillips/moc, krvi28/MOC_for_Rocket_Nozzle_Design

Design principle
----------------
The **wall is part of the characteristic mesh**. Each wall node inherits the
flow state of the previous interior node (wave cancellation) and is placed at
the intersection of:

1. the C+ characteristic from that interior node, and
2. the wall streamline from the previous wall node.

The full C+/C− net therefore fills the nozzle from throat to exit — the same
mesh used for contour design is drawn in MOC field plots.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from rockmoc.exceptions import DesignError, MeshError, ThermoError
from rockmoc.geometry.contour import Contour
from rockmoc.moc.equations import source_minus, source_plus
from rockmoc.moc.types import CharNet, FlowPoint, PointKind
from rockmoc.thermo.gas_protocol import MOCGas


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def _pm_state(nu: float, gas: MOCGas) -> tuple[float, float, float]:
    """(M, μ, ν_clamped) from Prandtl–Meyer angle."""
    nu_c = max(0.0, float(nu))
    if hasattr(gas, "prandtl_meyer_max"):
        nu_max = float(gas.prandtl_meyer_max()) * 0.999  # type: ignore[attr-defined]
    else:
        g = float(gas.gamma)
        gp = math.sqrt((g + 1.0) / (g - 1.0))
        nu_max = 0.5 * math.pi * (gp - 1.0) * 0.999
    if nu_c > nu_max:
        nu_c = nu_max
    if nu_c <= 1e-16:
        M = 1.01
        nu_c = float(gas.prandtl_meyer(M))
    else:
        try:
            M = float(gas.mach_from_prandtl_meyer(nu_c))
        except Exception:
            # Safe fallback if inverse fails near ν_max
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
    g_loc = float(getattr(gas, "gamma_at_mach", lambda _m: gas.gamma)(M))
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


def _n_points(n_lines: int) -> int:
    """Total characteristic points for N lines (centerline + wall inclusive)."""
    return int(n_lines + n_lines * 0.5 * (n_lines + 1))


def _mark_flags(n_lines: int, n_points: int) -> tuple[list[bool], list[bool]]:
    """Wall / centerline flags — literature index starts at 1."""
    on_wall = [False] * n_points
    on_cent = [False] * n_points
    j = 1 + n_lines
    k = 0
    for i in range(n_points):
        idx = i + 1
        if idx == j + k:
            on_wall[i] = True
            k += 1
            j += n_lines - k
    on_cent[0] = True
    for i in range(1, n_points):
        if on_wall[i - 1]:
            on_cent[i] = True
    return on_wall, on_cent


# ---------------------------------------------------------------------------
# Full MLN net (planar topology + optional axisym source correction)
# ---------------------------------------------------------------------------

def build_mln_full_net(
    gas: MOCGas,
    *,
    Rt: float,
    theta_max: float,
    n_chars: int,
    throat_mach: float = 1.01,
    delta: float = 1.0,
    transition: bool = True,
) -> tuple[CharNet, list[list[FlowPoint]]]:
    """Build full MLN characteristic net including wall nodes.

    Parameters
    ----------
    n_chars:
        Number of characteristic lines from the throat corner (N ≥ 4).
    delta:
        0 = planar (pure Riemann invariants), 1 = axisymmetric (source
        correction on invariants after each geometric step).
    """
    if n_chars < 4:
        raise DesignError("n_chars must be ≥ 4")
    if theta_max <= 0.0 or theta_max > math.radians(70.0):
        raise DesignError(f"theta_max out of range: {math.degrees(theta_max):.2f} deg")
    if Rt <= 0.0:
        raise DesignError("throat_radius must be positive")

    n = int(n_chars)
    # Equal flow-angle divisions 0 … θ_max (inclusive endpoints)
    dth = theta_max / (n - 1) if n > 1 else theta_max
    divs = [dth * i for i in range(n)]

    n_pts = _n_points(n)
    on_wall, on_cent = _mark_flags(n, n_pts)
    planar = abs(delta) < 0.5

    flow = [0.0] * n_pts
    pran = [0.0] * n_pts
    mach = [0.0] * n_pts
    mang = [0.0] * n_pts
    kneg = [0.0] * n_pts
    kpos = [0.0] * n_pts
    xx = [0.0] * n_pts
    yy = [0.0] * n_pts

    def _set_props(i: int, th: float, nu: float) -> None:
        M, mu, nu_c = _pm_state(nu, gas)
        flow[i] = th
        pran[i] = nu_c
        mach[i] = M
        mang[i] = mu
        kneg[i] = th + nu_c
        kpos[i] = th - nu_c

    def _axisym_interior(i: int, top: tuple[float, float, float, float, float], bot: int) -> None:
        """Axisymmetric interior unit process (averaged-coefficient predictor–corrector).

        ``top`` = (x, r, θ, ν, μ) of the upstream point on the C−; ``bot`` indexes
        the upstream point on the C+. Point ``i`` holds the planar predictor.
        """
        xt, rt, tht, nut, mut = top
        guess = (xx[i], yy[i], flow[i], pran[i])
        for _ in range(_UNIT_ITER):
            th_m, mu_m, r_m = 0.5 * (tht + flow[i]), 0.5 * (mut + mang[i]), 0.5 * (rt + yy[i])
            th_p, mu_p = 0.5 * (flow[bot] + flow[i]), 0.5 * (mang[bot] + mang[i])
            r_p = 0.5 * (yy[bot] + yy[i])
            x_new, y_new = _find_xy((xt, rt), (xx[bot], yy[bot]), th_m - mu_m, th_p + mu_p)
            kn = (tht + nut) + source_minus(th_m, mu_m, r_m, delta=delta) * (x_new - xt)
            kp = (flow[bot] - pran[bot]) + source_plus(th_p, mu_p, r_p, delta=delta) * (x_new - xx[bot])
            th_old, nu_old = flow[i], pran[i]
            xx[i], yy[i] = x_new, y_new
            th_new = 0.5 * (kn + kp)
            if not (math.isfinite(th_new) and math.isfinite(x_new) and abs(th_new) < 0.5 * math.pi):
                # Corrector diverged (near-sonic seed close to the axis): keep planar predictor
                xx[i], yy[i] = guess[0], guess[1]
                _set_props(i, guess[2], guess[3])
                break
            _set_props(i, th_new, max(0.5 * (kn - kp), 0.0))
            if abs(flow[i] - th_old) < _UNIT_TOL and abs(pran[i] - nu_old) < _UNIT_TOL:
                break

    def _axisym_axis(i: int, top: int, cnt: int) -> None:
        """Axis point: C− from ``top`` to r=0, θ=0; (θ+ν) carries the S− source."""
        for _ in range(_UNIT_ITER):
            th_m, mu_m, r_m = 0.5 * flow[top], 0.5 * (mang[top] + mang[i]), 0.5 * yy[top]
            x_new, _ = _find_xy((xx[top], yy[top]), (xx[cnt], 0.0), th_m - mu_m, 0.0)
            nu = kneg[top] + source_minus(th_m, mu_m, r_m, delta=delta) * (x_new - xx[top])
            nu_old = pran[i]
            xx[i] = max(x_new, xx[cnt] + 1e-12)
            _set_props(i, 0.0, max(nu, 0.0))
            if abs(pran[i] - nu_old) < _UNIT_TOL:
                break

    # ----- Point 0: first centerline node -----
    M0 = max(float(throat_mach), 1.01)
    nu0 = float(gas.prandtl_meyer(M0))
    mu0 = float(gas.mach_angle(M0))
    flow[0] = 0.0
    pran[0] = nu0 * 0.0  # classic MLN: start ν≈0 at near-sonic axis seed
    # Use M=1.01 so μ is defined; ν kept ~0 for Riemann start (nategphillips)
    mach[0] = M0
    mang[0] = mu0
    pran[0] = float(gas.prandtl_meyer(M0)) * 0.0
    # Recompute with tiny ν for numerical PM inverse stability
    _set_props(0, 0.0, 1e-8)
    mach[0] = M0
    mang[0] = mu0
    xx[0] = Rt / math.tan(mang[0] - flow[0])
    yy[0] = 0.0

    # ----- Axisym corner region: C+ rows launched from the initial line -----
    # The axis-reflected rows never reach the corner, and in axisymmetric flow
    # the corner region is not a simple wave, so cover it with C+ rows that
    # start on ray 0 (uniform near-sonic line) between point 0 and the corner.
    pre_rows: list[list[FlowPoint]] = []
    if not planar:
        from rockmoc.moc.mesh import interior_point

        nu_s = float(gas.prandtl_meyer(M0))
        m_pre = max(n // 2, 4)
        try:
            for p in range(m_pre):
                s_frac = ((p + 1) / (m_pre + 1)) ** 2  # clustered toward the corner
                cur = [_fp(s_frac * xx[0], Rt * (1.0 - s_frac), 0.0, nu_s, gas, PointKind.INTERIOR)]
                for j in range(1, n):
                    a = pre_rows[-1][j] if pre_rows else _fp(0.0, Rt, divs[j], divs[j], gas, PointKind.FAN)
                    cur.append(interior_point(a, cur[-1], gas, delta=delta))
                pre_rows.append(cur)
        except (MeshError, ThermoError):
            pre_rows = []

    # ----- Points 1 … N: first C−/C+ wave (ends at first wall) -----
    for i in range(1, n + 1):
        prv = i - 1
        if not on_wall[i]:
            th = divs[i]
            _set_props(i, th, th)  # simple-wave fan: ν = θ
            c_neg = flow[i] - mang[i]
            c_pos = 0.5 * ((flow[prv] + mang[prv]) + (flow[i] + mang[i]))
            xx[i], yy[i] = _find_xy((0.0, Rt), (xx[prv], yy[prv]), c_neg, c_pos)
            if not planar:
                # Corner state on this fan ray: θ = ν = divs[i] (PM valid at r = Rt)
                if pre_rows:
                    q = pre_rows[-1][i]
                    _axisym_interior(i, (q.x, q.r, q.theta, q.nu, q.mu), prv)
                else:
                    _, mu_c, _ = _pm_state(divs[i], gas)
                    _axisym_interior(i, (0.0, Rt, divs[i], divs[i], mu_c), prv)
        else:
            # First wall: inherit previous interior (cancellation)
            _set_props(i, flow[prv], pran[prv])
            c_neg = theta_max  # wall slope at corner region
            c_pos = 0.5 * ((flow[prv] + mang[prv]) + (flow[i] + mang[i]))
            xx[i], yy[i] = _find_xy((0.0, Rt), (xx[prv], yy[prv]), c_neg, c_pos)

    # ----- Remaining points -----
    jj = 0
    for i in range(n + 1, n_pts):
        prv = i - 1
        top = i - (n - jj)
        cnt = i - (n - jj) - 1

        if on_cent[i]:
            # Centerline: θ=0, K− from upper point
            kneg[i] = kneg[top]
            th = 0.0
            nu = kneg[i] - th
            _set_props(i, th, max(nu, 0.0))
            c_neg = 0.5 * ((flow[top] - mang[top]) + (flow[i] - mang[i]))
            c_pos = 0.0
            xx[i], yy[i] = _find_xy((xx[top], yy[top]), (xx[cnt], yy[cnt]), c_neg, c_pos)
            yy[i] = 0.0
            if xx[i] < xx[cnt]:
                xx[i] = xx[cnt] + 1e-12
            if not planar:
                _axisym_axis(i, top, cnt)

        elif not on_wall[i]:
            # Interior: K− from top, K+ from previous
            kneg[i] = kneg[top]
            kpos[i] = kpos[prv]
            th = 0.5 * (kneg[i] + kpos[i])
            nu = 0.5 * (kneg[i] - kpos[i])
            _set_props(i, th, max(nu, 0.0))
            c_neg = 0.5 * ((flow[top] - mang[top]) + (flow[i] - mang[i]))
            c_pos = 0.5 * ((flow[prv] + mang[prv]) + (flow[i] + mang[i]))
            xx[i], yy[i] = _find_xy((xx[top], yy[top]), (xx[prv], yy[prv]), c_neg, c_pos)
            if not planar:
                _axisym_interior(i, (xx[top], yy[top], flow[top], pran[top], mang[top]), prv)
            if yy[i] < 0.0:
                yy[i] = 0.0

        else:
            # Wall: cancel — inherit previous interior state
            _set_props(i, flow[prv], pran[prv])
            # Streamline slope = average wall flow angles
            c_neg = 0.5 * (flow[top] + flow[i])  # top is previous wall
            c_pos = 0.5 * ((flow[prv] + mang[prv]) + (flow[i] + mang[i]))
            xx[i], yy[i] = _find_xy((xx[top], yy[top]), (xx[prv], yy[prv]), c_neg, c_pos)
            if yy[i] < yy[top]:
                # Force mild expansion if fold
                yy[i] = yy[top] + 1e-9 * Rt
                if abs(math.tan(c_neg)) > 1e-12:
                    xx[i] = xx[top] + (yy[i] - yy[top]) / math.tan(c_neg)
            if xx[i] < xx[top]:
                xx[i] = xx[top] + 1e-12
            jj += 1

    # ----- Pack into CharNet + structured rows -----
    net = CharNet()
    throat = _fp(0.0, Rt, 0.0, 1e-8, gas, PointKind.THROAT)
    net.add(throat)

    # Fan property seeds (corner)
    net.fan = [_fp(0.0, Rt, divs[j], divs[j], gas, PointKind.FAN) for j in range(n)]

    # Build rows: each wave from centerline to wall
    # Reconstruct by walking points in order
    rows: list[list[FlowPoint]] = []
    axis_list: list[FlowPoint] = []
    wall_list: list[FlowPoint] = [throat]
    current_row: list[FlowPoint] = []

    for i in range(n_pts):
        if on_cent[i]:
            if current_row:
                rows.append(current_row)
            kind = PointKind.AXIS
            p = _fp(xx[i], 0.0, flow[i], pran[i], gas, kind)
            net.add(p)
            axis_list.append(p)
            current_row = [p]
        elif on_wall[i]:
            kind = PointKind.WALL if i < n_pts - 1 else PointKind.EXIT
            p = _fp(xx[i], yy[i], flow[i], pran[i], gas, kind)
            net.add(p)
            wall_list.append(p)
            current_row.append(p)
            rows.append(current_row)
            current_row = []
        else:
            p = _fp(xx[i], yy[i], flow[i], pran[i], gas, PointKind.INTERIOR)
            net.add(p)
            current_row.append(p)

    if current_row:
        rows.append(current_row)

    # Clean wall angles: non-increasing magnitude toward axial exit
    if wall_list:
        for wi in range(1, len(wall_list)):
            th_w = wall_list[wi].theta
            # Clamp small negatives from numerical noise
            if th_w < 0.0:
                wall_list[wi] = wall_list[wi].copy(theta=0.0)
        w = wall_list[-1]
        wall_list[-1] = _fp(w.x, w.r, 0.0, w.nu, gas, PointKind.EXIT)
        if rows and rows[-1]:
            rows[-1][-1] = wall_list[-1]

    net.axis = axis_list
    net.wall = wall_list
    net.interior = rows
    if not planar and transition:
        _axisym_transition(gas, net, Rt=Rt, theta_max=theta_max, delta=delta, pre_rows=pre_rows)
        rows = net.interior
    return net, rows


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


def _transition_wall(
    gas: MOCGas,
    bases: list[list[FlowPoint]],
    *,
    Rt: float,
    delta: float,
    mdot: float | None = None,
    start: FlowPoint | None = None,
) -> tuple[list[FlowPoint], list[list[FlowPoint]]]:
    """Transition region + mass-flow wall (Zucrow & Hoffman).

    ``bases`` are C+ polylines ordered from the axis point E (where M = Me)
    outward; each ends on the last expansion C−. The transition region is the
    Goursat problem between that C− and the uniform exit C+ from E (M = Me,
    θ = 0). The wall is the streamline carrying the sonic-throat mass flow,
    located along each C+ (``mdot`` defaults to the uniform sonic throat).
    Returns (wall from upstream to exit, rows).
    """
    from rockmoc.moc.mesh import interior_point

    planar = abs(delta) < 0.5
    E = bases[0][0]
    if planar and start is not None:
        # Planar: downstream of the last C− the flow is a simple wave, so each
        # C+ is straight with the state of its base point — exact streamline march.
        wall_p: list[FlowPoint] = []
        rows_p: list[list[FlowPoint]] = []
        prev_w = start
        for row in reversed(bases):
            b = row[-1]
            th_w = 0.0 if b is E else b.theta
            x, y = _find_xy((prev_w.x, prev_w.r), (b.x, b.r), 0.5 * (prev_w.theta + th_w), b.theta + b.mu)
            w = _fp(x, y, th_w, b.nu, gas, PointKind.EXIT if b is E else PointKind.WALL)
            wall_p.append(w)
            rows_p.append(row + [w])
            prev_w = w
        return wall_p, rows_p[::-1]
    n_u = max(len(bases), 8)
    if mdot is None:
        mdot = Rt if planar else math.pi * Rt * Rt
    ar_e = float(gas.area_ratio(E.M))
    r_e = mdot * ar_e if planar else math.sqrt(mdot * ar_e / math.pi)
    t_e = r_e / math.sin(E.mu)
    prev = [
        FlowPoint(x=E.x + t * math.cos(E.mu), r=t * math.sin(E.mu), theta=0.0, nu=E.nu,
                  M=E.M, mu=E.mu, kind=PointKind.INTERIOR, gamma_local=E.gamma_local)
        for t in np.linspace(0.0, t_e, n_u + 1)
    ]
    ext_rows = [prev]
    for row in bases[1:]:
        cur = [row[-1]]
        for i in range(1, n_u + 1):
            cur.append(interior_point(prev[i], cur[-1], gas, delta=delta))
        ext_rows.append(cur)
        prev = cur

    wall: list[FlowPoint] = []
    rows: list[list[FlowPoint]] = []
    for row, ext in zip(bases, ext_rows):
        line = row + ext[1:]
        m = _poly_mass(line, gas, planar=planar)
        j = next((j for j in range(1, len(m)) if m[j] >= mdot), len(m) - 1)
        f = min(max((mdot - m[j - 1]) / max(m[j] - m[j - 1], 1e-300), 0.0), 1.0)
        a, b = line[j - 1], line[j]
        w = _fp(a.x + f * (b.x - a.x), a.r + f * (b.r - a.r), a.theta + f * (b.theta - a.theta),
                a.nu + f * (b.nu - a.nu), gas, PointKind.WALL)
        wall.append(w)
        rows.append(line[:j] + [w])
    wall.reverse()
    rows.reverse()
    wall[-1] = _fp(wall[-1].x, r_e, 0.0, E.nu, gas, PointKind.EXIT)
    rows[-1][-1] = wall[-1]
    return wall, rows


def _axisym_transition(
    gas: MOCGas,
    net: CharNet,
    *,
    Rt: float,
    theta_max: float,
    delta: float,
    pre_rows: list[list[FlowPoint]] | None = None,
) -> None:
    """Sharp-corner MLN: transition region + mass-flow wall from the kernel rows."""
    kernel = [row[:-1] for row in net.interior]  # drop provisional wall nodes
    # Corner-region rows prepend axis point 0 → start so the mass integral
    # covers the uniform strip of the initial line.
    bases = kernel[::-1] + [[kernel[0][0]] + row for row in (pre_rows or [])[::-1]]
    wall, rows = _transition_wall(gas, bases, Rt=Rt, delta=delta)
    throat = _fp(0.0, Rt, theta_max, theta_max, gas, PointKind.THROAT)
    net.points = [throat] + list(net.fan) + [p for row in rows for p in row]
    net.wall = [throat] + wall
    net.interior = rows


def build_axisym_kernel(
    gas: MOCGas,
    *,
    Rt: float,
    theta_max: float,
    n_chars: int,
    throat_mach: float = 1.01,
    delta: float = 1.0,
    unit_tol: float = 1e-9,
    unit_iter: int = 12,
) -> tuple[CharNet, list[list[FlowPoint | None]], list[FlowPoint]]:
    """Compatibility wrapper returning ``(net, grid, fan)``."""
    _ = unit_tol, unit_iter
    net, rows = build_mln_full_net(
        gas,
        Rt=Rt,
        theta_max=theta_max,
        n_chars=n_chars,
        throat_mach=throat_mach,
        delta=delta,
    )
    grid: list[list[FlowPoint | None]] = [list(row) for row in rows]  # type: ignore[list-item]
    return net, grid, list(net.fan)


# ---------------------------------------------------------------------------
# Smooth-throat (circular arc) ideal contour net
# ---------------------------------------------------------------------------

def build_ideal_net(
    gas: MOCGas,
    *,
    Rt: float,
    Rc: float,
    nu_target: float,
    n_chars: int,
    throat_mach: float = 1.01,
    delta: float = 1.0,
    transition: bool = True,
    initial_line: Literal["sauer", "uniform"] = "sauer",
    attach_offset: int = 0,
) -> CharNet:
    """Ideal (uniform-exit) contour with a circular-arc throat of radius Rc.

    ``initial_line="sauer"`` starts from Sauer's transonic solution on the
    vertical line where the axis reaches throat_mach; ``"uniform"`` starts
    from uniform θ = 0, M = throat_mach on the throat plane (cruder, but lets
    the arc end before the Sauer line for small Rc / low Me).

    Characteristic lattice P(α, β): C+ lines α, C− lines β, marched along
    diagonals s = α + β from the initial line x = 0 (uniform θ = 0,
    M = throat_mach), whose points are I_j = P(m − j, j).

    * interior P(α, β) from C− through P(α−1, β) and C+ through P(α, β−1)
    * wall    P(α, α+m) on the arc x = Rc sinφ, r = Rt + Rc(1 − cosφ):
      C+ ∩ circle, θ_W = φ, ν_W from C+ compatibility
    * axis    P(β+m, β) from the C− through P(β+m−1, β)

    The lattice upstream of any C− does not depend on where the arc ends, so
    the arc is marched until the C− from a wall node reaches the axis at
    ν ≥ ``nu_target``; that node is the attachment (θ_a = its angle) and E is
    the axis end of its C−. ``attach_offset = 1`` attaches one wall node
    earlier (E just below the target) so callers can interpolate between the
    two. The transition region + mass-flow wall follow (:func:`_transition_wall`).
    """
    from rockmoc.moc.mesh import axis_point, interior_point

    m = max(int(n_chars), 6)
    yc = Rt + Rc
    # Sauer transonic solution (Zucrow & Hoffman) for the initial data line:
    # the vertical line through the axial point where M = throat_mach. A
    # vertical line stays spacelike near M = 1 (iso-Mach lines do not).
    g = float(gas.gamma)
    dl = 0.0 if abs(delta) < 0.5 else 1.0
    alpha = math.sqrt((1.0 + dl) / ((g + 1.0) * Rc * Rt))
    k2 = (g + 1.0) * alpha**2 / (2.0 * (1.0 + dl))
    x_shift = (g + 1.0) * alpha * Rt * Rt / (2.0 * (3.0 + dl))  # axial sonic point − throat plane
    M0 = max(float(throat_mach), 1.0 + 1e-6)
    xs = (math.sqrt(0.5 * (g + 1.0) * M0 * M0 / (1.0 + 0.5 * (g - 1.0) * M0 * M0)) - 1.0) / alpha
    x_line = xs + x_shift if initial_line == "sauer" else 0.0
    if x_line >= Rc:
        raise MeshError("throat_mach too high for this throat radius of curvature")
    phi0 = math.asin(x_line / Rc)
    r_w = Rt + Rc * (1.0 - math.cos(phi0))
    m_star_max = 0.999 * math.sqrt((g + 1.0) / (g - 1.0))
    P: dict[tuple[int, int], FlowPoint] = {}
    for j in range(m + 1):
        r = r_w * j / m
        u = alpha * xs + k2 * r * r
        v = ((g + 1.0) * alpha**2 * xs * r / (1.0 + dl)
             + (g + 1.0) ** 2 * alpha**3 * r**3 / (2.0 * (1.0 + dl) * (3.0 + dl)))
        ms = min(1.0 + u, m_star_max)
        M = math.sqrt(2.0 * ms * ms / ((g + 1.0) - (g - 1.0) * ms * ms))
        th = phi0 if j == m else math.atan2(v, 1.0 + u)
        if initial_line != "sauer":
            M, th = M0, 0.0
        kind = PointKind.AXIS if j == 0 else (PointKind.THROAT if j == m else PointKind.INTERIOR)
        P[(m - j, j)] = _fp(x_line, r, th, float(gas.prandtl_meyer(M)), gas, kind)

    def _arc_wall(b: FlowPoint) -> FlowPoint:
        al = b.theta + b.mu
        w = None
        for _ in range(_UNIT_ITER):
            dx, dr = math.cos(al), math.sin(al)
            ox, orr = b.x, b.r - yc
            bq = 2.0 * (ox * dx + orr * dr)
            disc = bq * bq - 4.0 * (ox * ox + orr * orr - Rc * Rc)
            if disc < 0.0:
                raise MeshError("C+ misses the throat arc")
            t = 0.5 * (-bq - math.sqrt(disc))
            x, r = b.x + t * dx, b.r + t * dr
            phi = math.atan2(x, yc - r)
            th_p = 0.5 * (b.theta + phi)
            mu_p = b.mu if w is None else 0.5 * (b.mu + w.mu)
            kp = (b.theta - b.nu) + source_plus(th_p, mu_p, 0.5 * (b.r + r), delta=delta) * (x - b.x)
            w_new = _fp(x, r, phi, max(phi - kp, 0.0), gas, PointKind.WALL)
            if w is not None and abs(w_new.nu - w.nu) < _UNIT_TOL:
                return w_new
            w = w_new
            al = 0.5 * (b.theta + b.mu + w.theta + w.mu)
        return w  # type: ignore[return-value]

    def _point(a: int, b: int) -> FlowPoint:
        if b - a == m:
            return _arc_wall(P[(a, b - 1)])
        if a - b == m:
            return axis_point(P[(a - 1, b)], gas, delta=delta)
        return interior_point(P[(a - 1, b)], P[(a, b - 1)], gas, delta=delta)

    # March diagonals s = α + β with the arc continuing; stop when the C− from
    # a wall node (β = b, wall at α = b − m) reaches the axis at ν ≥ ν_target.
    b_end = None
    for s_ in range(m + 1, m + 1 + 400 * m):
        for a in range(max((s_ - m + 1) // 2, 0), (s_ + m) // 2 + 1):
            b = s_ - a
            if abs(a - b) > m or (a, b) in P:
                continue
            try:
                P[(a, b)] = _point(a, b)
            except (MeshError, ThermoError, KeyError):
                continue  # far ahead on the arc; only needed if E is not reached first
            if a - b == m and b >= m and P[(a, b)].nu >= nu_target:
                b_end = b
                break
        if b_end is not None:
            break
        if not any(k[0] + k[1] == s_ for k in P):
            raise MeshError("characteristic net stopped before reaching the target Mach")
    if b_end is None:
        raise MeshError("target Mach not reached on the axis")
    b_end -= attach_offset
    if b_end <= m:
        raise MeshError("attachment upstream of the initial data line")

    a_E = b_end + m
    a_star = b_end - m
    P = {k: v for k, v in P.items() if k[1] <= b_end and k[0] <= a_E}
    wall = [P[(aw, aw + m)] for aw in range(0, a_star + 1)]

    net = CharNet()
    rows: list[list[FlowPoint]] = []
    for a in sorted({k[0] for k in P}):
        rows.append([P[k] for k in sorted((k for k in P if k[0] == a), key=lambda k: k[1])])
    net.axis = [P[(bb + m, bb)] for bb in range(0, b_end + 1) if (bb + m, bb) in P]
    net.wall = list(wall)
    net.interior = rows
    if transition:
        bases: list[list[FlowPoint]] = []
        initial = [P[(m - j, j)] for j in range(m + 1)]
        for a in range(a_E, a_star, -1):
            line = [P[(a, bb)] for bb in range(0, b_end + 1) if (a, bb) in P]
            if a <= m:  # C+ starts on the initial line: include the strip below it
                line = initial[: m - a] + line
            bases.append(line)
        # Mass across the last C− (E → attachment): the wall streamline then
        # passes exactly through the attachment node.
        c_minus = [P[(a, b_end)] for a in range(a_E, a_star - 1, -1)]
        mdot = _poly_mass(c_minus, gas, planar=abs(delta) < 0.5)[-1]
        t_wall, t_rows = _transition_wall(gas, bases, Rt=Rt, delta=delta, mdot=mdot, start=wall[-1])
        net.wall = list(wall) + [w for w in t_wall if w.x > wall[-1].x]
        net.interior = rows + t_rows
    net.points = [p for r_ in net.interior for p in r_]
    return net


# ---------------------------------------------------------------------------
# Contour helpers
# ---------------------------------------------------------------------------

def _wall_to_contour(wall: list[FlowPoint], Rt: float) -> Contour:
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


# ---------------------------------------------------------------------------
# Public design API
# ---------------------------------------------------------------------------

@dataclass
class AxisymHFResult:
    kind: Literal["mln", "ideal", "tic"]
    exit_mach_design: float
    gamma: float
    throat_radius: float
    n_chars: int
    theta_max_rad: float
    contour: Contour
    net: CharNet | None
    area_ratio_geometric: float
    area_ratio_1d: float
    length: float
    max_wall_angle_rad: float
    axis_exit_mach: float | None
    meta: dict = field(default_factory=dict)

    def summary(self) -> dict:
        eps_g = self.area_ratio_geometric
        eps_1 = self.area_ratio_1d
        return {
            "method": f"axisym_hf_{self.kind}",
            "exit_mach_design": self.exit_mach_design,
            "gamma": self.gamma,
            "throat_radius_m": self.throat_radius,
            "n_chars": self.n_chars,
            "theta_max_deg": math.degrees(self.theta_max_rad),
            "area_ratio_geometric": eps_g,
            "area_ratio_1d": eps_1,
            "area_ratio_rel_error": abs(eps_g - eps_1) / eps_1 if eps_1 > 0 else None,
            "length_m": self.length,
            "length_over_Rt": self.length / self.throat_radius,
            "max_wall_angle_deg": math.degrees(self.max_wall_angle_rad),
            "axis_exit_mach": self.axis_exit_mach,
            "exit_radius_m": float(self.contour.r[-1]),
            **self.meta,
        }


def design_axisym_mln_hf(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    n_chars: int = 40,
    throat_mach: float = 1.01,
    delta: float = 1.0,
    match_exit_area: bool = True,
    store_net: bool = True,
    resample: int | None = None,
) -> AxisymHFResult:
    """Minimum-length nozzle via full MOC net (planar or axisymmetric)."""
    if exit_mach <= 1.0:
        raise DesignError("exit_mach must be > 1")
    if throat_radius <= 0.0:
        raise DesignError("throat_radius must be positive")
    if n_chars < 4:
        raise DesignError("n_chars must be ≥ 4")

    Rt = float(throat_radius)
    Me = float(exit_mach)
    nu_e = float(gas.prandtl_meyer(Me))
    theta_max = 0.5 * nu_e

    def _net(th: float, transition: bool = True) -> CharNet:
        return build_mln_full_net(
            gas, Rt=Rt, theta_max=th, n_chars=n_chars, throat_mach=throat_mach, delta=delta,
            transition=transition,
        )[0]

    if abs(delta) > 0.5:
        # Axisymmetric: sources raise axis ν above 2θ_max, so solve for the
        # corner angle that puts the kernel's last axis point exactly at ν_e.
        from scipy.optimize import brentq

        def _resid(th: float) -> float:
            try:
                return _net(th, transition=False).axis[-1].nu - nu_e
            except (MeshError, ThermoError):
                return nu_e  # corner too steep: over-expands past ν_max

        lo, hi = 0.2 * nu_e, 0.5 * nu_e
        if _resid(lo) < 0.0 < _resid(hi):
            theta_max = float(brentq(_resid, lo, hi, xtol=1e-10))

    net = _net(theta_max)

    # Wall angles should start near θ_max: set throat wall node angle
    if net.wall and len(net.wall) >= 2:
        net.wall[0] = _fp(0.0, Rt, theta_max, theta_max, gas, PointKind.THROAT)

    contour = _wall_to_contour(net.wall, Rt)
    # Force first wall angle to θ_max for contour colouring
    if contour.theta is not None and len(contour.theta):
        th = contour.theta.copy()
        th[0] = theta_max
        th[-1] = 0.0
        contour = Contour(x=contour.x, r=contour.r, theta=th, M=contour.M)

    eps_1d = float(gas.area_ratio(Me))
    if abs(delta) < 0.5:
        r_target = Rt * eps_1d
    else:
        r_target = Rt * math.sqrt(eps_1d)

    if match_exit_area and contour.r[-1] > Rt + 1e-15:
        scale = (r_target - Rt) / (float(contour.r[-1]) - Rt)
        contour = _isotropic_match_exit(contour, Rt, r_target)
        _scale_net(net, Rt, scale, x0=0.0)
        if net.wall:
            net.wall[0] = _fp(0.0, Rt, theta_max, theta_max, gas, PointKind.THROAT)
            last = net.wall[-1]
            net.wall[-1] = _fp(last.x, r_target, 0.0, last.nu, gas, PointKind.EXIT)

    n_out = resample if resample is not None else max(80, n_chars * 3)
    contour_out = contour.resample(n_out)

    if abs(delta) < 0.5:
        eps_g = float(contour_out.r[-1] / Rt)
    else:
        eps_g = float((contour_out.r[-1] / Rt) ** 2)

    axis_M = net.axis[-1].M if net.axis else None
    max_ang = float(np.max(np.abs(contour_out.theta))) if contour_out.theta is not None else theta_max

    return AxisymHFResult(
        kind="mln",
        exit_mach_design=Me,
        gamma=float(gas.gamma),
        throat_radius=Rt,
        n_chars=n_chars,
        theta_max_rad=theta_max,
        contour=contour_out,
        net=net if store_net else None,
        area_ratio_geometric=eps_g,
        area_ratio_1d=eps_1d,
        length=float(contour_out.length),
        max_wall_angle_rad=max_ang,
        axis_exit_mach=float(axis_M) if axis_M is not None else None,
        meta={
            "delta": delta,
            "match_exit_area": match_exit_area,
            "nu_e_rad": nu_e,
            "throat_mach": throat_mach,
            "n_wall_nodes": len(net.wall),
            "n_net_points": len(net.points),
            "algorithm": "full_moc_net_mln",
        },
    )


def _attach_angle(net: CharNet, Rt: float, Rc: float) -> float:
    """Angle of the last wall node lying on the throat arc."""
    on_arc = [p for p in net.wall if abs(math.hypot(p.x, Rt + Rc - p.r) - Rc) < 1e-9 * Rc]
    return on_arc[-1].theta if on_arc else 0.0


def design_axisym_ideal_hf(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    n_chars: int = 40,
    throat_rc_over_rt: float = 1.5,
    throat_mach: float = 1.01,
    delta: float = 1.0,
    match_exit_area: bool = True,
    store_net: bool = True,
) -> AxisymHFResult:
    """Ideal (uniform, axial exit) contour with a circular-arc throat — full MOC.

    The arc attachment angle θ_a is solved so that the last expansion C−
    reaches the axis at M = Me; the downstream wall closes the throat mass flow.
    """
    if exit_mach <= 1.0:
        raise DesignError("exit_mach must be > 1")
    if throat_radius <= 0.0 or throat_rc_over_rt <= 0.0:
        raise DesignError("throat_radius and throat_rc_over_rt must be positive")
    Rt = float(throat_radius)
    Rc = throat_rc_over_rt * Rt
    Me = float(exit_mach)
    nu_e = float(gas.prandtl_meyer(Me))

    net = net_lo = None
    for start in ("sauer", "uniform"):
        try:
            kw = dict(gas=gas, Rt=Rt, Rc=Rc, nu_target=nu_e, n_chars=n_chars,
                      throat_mach=throat_mach, delta=delta, initial_line=start)
            net = build_ideal_net(**kw)  # type: ignore[arg-type]
            net_lo = build_ideal_net(**kw, attach_offset=1)  # type: ignore[arg-type]
            break
        except (MeshError, ThermoError):
            net = net_lo = None
    if net is None or net_lo is None:
        raise DesignError(
            f"ideal contour not found for Me={Me}; "
            "try a larger throat_rc_over_rt or lower throat_mach"
        )
    theta_a = _attach_angle(net, Rt, Rc)

    # The attachment is discrete: blend the two bracketing designs linearly in ν_E
    nu_hi, nu_lo = net.axis[-1].nu, net_lo.axis[-1].nu
    wgt = min(max((nu_e - nu_lo) / max(nu_hi - nu_lo, 1e-12), 0.0), 1.0)
    n_bl = max(200, 6 * n_chars)
    c_hi = _wall_to_contour(net.wall, Rt).resample(n_bl)
    c_lo = _wall_to_contour(net_lo.wall, Rt).resample(n_bl)
    th = (1.0 - wgt) * c_lo.theta + wgt * c_hi.theta
    th[0] = 0.0
    th[-1] = 0.0
    contour = Contour(
        x=(1.0 - wgt) * c_lo.x + wgt * c_hi.x,
        r=(1.0 - wgt) * c_lo.r + wgt * c_hi.r,
        theta=th,
        M=(1.0 - wgt) * c_lo.M + wgt * c_hi.M,
    )
    theta_a = (1.0 - wgt) * _attach_angle(net_lo, Rt, Rc) + wgt * theta_a
    eps_1d = float(gas.area_ratio(Me))
    planar = abs(delta) < 0.5
    r_t = Rt * (eps_1d if planar else math.sqrt(eps_1d))
    if match_exit_area and contour.r[-1] > Rt:
        contour = _isotropic_match_exit(contour, Rt, r_t)
    contour = contour.resample(max(100, n_chars * 3))
    ratio = float(contour.r[-1] / Rt)
    eps_g = ratio if planar else ratio**2
    max_ang = float(np.max(np.abs(contour.theta)))

    return AxisymHFResult(
        kind="ideal",
        exit_mach_design=Me,
        gamma=float(gas.gamma),
        throat_radius=Rt,
        n_chars=n_chars,
        theta_max_rad=theta_a,
        contour=contour,
        net=net if store_net else None,
        area_ratio_geometric=eps_g,
        area_ratio_1d=eps_1d,
        length=float(contour.length),
        max_wall_angle_rad=max_ang,
        axis_exit_mach=float(net.axis[-1].M),
        meta={
            "delta": delta,
            "nu_e_rad": nu_e,
            "throat_mach": throat_mach,
            "throat_rc_over_rt": throat_rc_over_rt,
            "theta_attach_deg": math.degrees(theta_a),
            "initial_line": start,
            "match_exit_area": match_exit_area,
            "n_wall_nodes": len(net.wall),
            "algorithm": "full_moc_net_ideal_arc",
        },
    )


def design_axisym_tic_hf(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    truncate_fraction: float = 0.80,
    n_chars: int = 40,
    **kwargs,
) -> AxisymHFResult:
    """Truncated ideal contour at a fraction of ideal length."""
    if not (0.05 < truncate_fraction <= 1.0):
        raise DesignError("truncate_fraction must be in (0.05, 1]")
    store = kwargs.pop("store_net", True)
    full = design_axisym_ideal_hf(
        gas,
        throat_radius=throat_radius,
        exit_mach=exit_mach,
        n_chars=n_chars,
        store_net=store,
        **kwargs,
    )
    contour = full.contour.truncate(truncate_fraction)
    ratio = float(contour.r[-1] / throat_radius)
    eps_g = ratio if abs(full.meta.get("delta", 1.0)) < 0.5 else ratio**2

    net = full.net
    if net is not None and store:
        x_cut = float(contour.x[-1])
        net.wall = [p for p in net.wall if p.x <= x_cut + 1e-12]
        if net.wall:
            th_e = float(contour.theta[-1]) if contour.theta is not None else net.wall[-1].theta
            Me_e = float(contour.M[-1]) if contour.M is not None else net.wall[-1].M
            nu_e = float(gas.prandtl_meyer(max(Me_e, 1.01)))
            net.wall[-1] = _fp(x_cut, float(contour.r[-1]), th_e, nu_e, gas, PointKind.EXIT)
        net.axis = [p for p in net.axis if p.x <= x_cut + 1e-9]
        net.interior = [[p for p in row if p.x <= x_cut + 1e-9] for row in net.interior]
        net.interior = [row for row in net.interior if len(row) >= 2]
        net.points = [p for p in net.points if p.x <= x_cut + 1e-9]

    max_ang = float(np.max(np.abs(contour.theta))) if contour.theta is not None else full.max_wall_angle_rad
    return AxisymHFResult(
        kind="tic",
        exit_mach_design=full.exit_mach_design,
        gamma=full.gamma,
        throat_radius=throat_radius,
        n_chars=n_chars,
        theta_max_rad=full.theta_max_rad,
        contour=contour,
        net=net,
        area_ratio_geometric=eps_g,
        area_ratio_1d=full.area_ratio_1d,
        length=float(contour.length),
        max_wall_angle_rad=max_ang,
        axis_exit_mach=full.axis_exit_mach,
        meta={
            **full.meta,
            "truncate_fraction": truncate_fraction,
            "ideal_length_m": full.length,
        },
    )


def mesh_convergence_axisym(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    n_list: list[int] | None = None,
    delta: float = 1.0,
    match_exit_area: bool = True,
) -> list[dict]:
    if n_list is None:
        n_list = [12, 20, 40, 60, 80]
    rows: list[dict] = []
    prev = None
    for n_c in n_list:
        r = design_axisym_mln_hf(
            gas,
            throat_radius=throat_radius,
            exit_mach=exit_mach,
            n_chars=n_c,
            delta=delta,
            match_exit_area=match_exit_area,
            store_net=False,
        )
        s = r.summary()
        row = {
            "n_chars": n_c,
            "length_over_Rt": s["length_over_Rt"],
            "area_ratio_geometric": s["area_ratio_geometric"],
            "area_ratio_1d": s["area_ratio_1d"],
            "area_ratio_rel_error": s["area_ratio_rel_error"],
            "theta_max_deg": s["theta_max_deg"],
            "axis_exit_mach": s["axis_exit_mach"],
        }
        if prev is not None:
            c0 = prev.contour.resample(100)
            c1 = r.contour.resample(100)
            r1i = np.interp(c0.x, c1.x, c1.r, left=c1.r[0], right=c1.r[-1])
            row["contour_L2_vs_prev"] = float(
                np.sqrt(np.mean((r1i - c0.r) ** 2)) / throat_radius
            )
        rows.append(row)
        prev = r
    return rows
