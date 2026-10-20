"""Sharp-corner minimum-length nozzle kernel (planar topology + axisymmetric unit processes)."""

from __future__ import annotations

import math

from rockmoc.exceptions import DesignError, MeshError, ThermoError, warn_fallback
from rockmoc.moc._common import _UNIT_ITER, _UNIT_TOL, _find_xy, _fp, _pm_state
from rockmoc.moc.equations import source_minus, source_plus
from rockmoc.moc.transition import _transition_wall
from rockmoc.moc.types import CharNet, FlowPoint, PointKind
from rockmoc.thermo.gas_protocol import MOCGas


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
                warn_fallback("axisymmetric corrector diverged near the axis; planar predictor kept")
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

    # ----- Point 0: first centerline node (end of fan ray 0) -----
    # M0 > 1 keeps μ0 < 90° so ray 0 reaches the axis at a finite x.
    M0 = max(float(throat_mach), 1.01)
    mu0 = float(gas.mach_angle(M0))
    # Axisymmetric: start the fan from a consistent M0 state (ν0 = ν(M0), fan
    # ray j carries ν0 + θ_j). Planar keeps the classic sonic start (ν0 = 0),
    # pinned by the NASA TM-X-1502 validation.
    nu0 = 0.0 if planar else float(gas.prandtl_meyer(M0))
    _set_props(0, 0.0, max(nu0, 1e-8))
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
                    a = pre_rows[-1][j] if pre_rows else _fp(0.0, Rt, divs[j], nu0 + divs[j], gas, PointKind.FAN)
                    cur.append(interior_point(a, cur[-1], gas, delta=delta))
                pre_rows.append(cur)
        except (MeshError, ThermoError):
            warn_fallback("corner-region characteristics failed; corner region left unresolved")
            pre_rows = []

    # ----- Points 1 … N: first C−/C+ wave (ends at first wall) -----
    for i in range(1, n + 1):
        prv = i - 1
        if not on_wall[i]:
            th = divs[i]
            _set_props(i, th, nu0 + th)  # simple-wave fan: ν = ν0 + θ
            c_neg = flow[i] - mang[i]
            c_pos = 0.5 * ((flow[prv] + mang[prv]) + (flow[i] + mang[i]))
            xx[i], yy[i] = _find_xy((0.0, Rt), (xx[prv], yy[prv]), c_neg, c_pos)
            if not planar:
                # Corner state on this fan ray: θ = ν = divs[i] (PM valid at r = Rt)
                if pre_rows:
                    q = pre_rows[-1][i]
                    _axisym_interior(i, (q.x, q.r, q.theta, q.nu, q.mu), prv)
                else:
                    _, mu_c, _ = _pm_state(nu0 + divs[i], gas)
                    _axisym_interior(i, (0.0, Rt, divs[i], nu0 + divs[i], mu_c), prv)
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
    net.fan = [_fp(0.0, Rt, divs[j], nu0 + divs[j], gas, PointKind.FAN) for j in range(n)]

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
    # C− lines = fan rays: corner → (corner-region rows) → kernel rows; ray m
    # meets kernel row k at index m − k (index 0 is the axis point).
    for m_ray in range(1, n):
        line = [throat] + [pr[m_ray] for pr in pre_rows]
        line += [rows[k][m_ray - k] for k in range(m_ray + 1)
                 if k < len(rows) and m_ray - k < len(rows[k]) - 1]
        net.c_minus.append(line)
    if not planar and transition:
        _axisym_transition(gas, net, Rt=Rt, theta_max=theta_max, delta=delta, pre_rows=pre_rows)
        rows = net.interior
    return net, rows


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
    wall, rows = _transition_wall(gas, bases, Rt=Rt, delta=delta, c_minus_out=net.c_minus)
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
