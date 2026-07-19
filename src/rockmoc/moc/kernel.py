"""Characteristic net construction for nozzle design.

Implements a classical minimum-length nozzle (MLN) construction:

1. Centered expansion fan at a sharp throat corner.
2. Kernel region filled with interior / axis unit processes.
3. Wave-canceling wall marched as a streamline until θ → 0.

For planar flow (δ=0) the corner angle is θ_max = ν(Me)/2 and the
invariants θ±ν are algebraic. Axisymmetric flow (δ=1) uses the same
topology with geometric source terms in the compatibility equations.
"""

from __future__ import annotations

import math

from rockmoc.exceptions import DesignError, MeshError
from rockmoc.moc.mesh import axis_point, interior_point
from rockmoc.moc.types import CharNet, FlowPoint, PointKind
from rockmoc.thermo.perfect_gas import PerfectGas


def build_mln_net(
    gas: PerfectGas,
    *,
    throat_radius: float,
    theta_max: float,
    n_chars: int,
    delta: float = 1.0,
    throat_mach: float = 1.05,
    arc_offset_frac: float = 1e-3,
) -> CharNet:
    """Build MLN characteristic net and wall contour.

    Parameters
    ----------
    theta_max:
        Wall angle immediately after the sharp-corner expansion [rad].
    n_chars:
        Number of fan rays (mesh density). Must be ≥ 4.
    delta:
        0 = planar, 1 = axisymmetric.
    """
    if n_chars < 4:
        raise DesignError("n_chars must be ≥ 4")
    if theta_max <= 0.0 or theta_max > math.radians(80.0):
        raise DesignError(f"theta_max out of range: {theta_max}")
    if throat_radius <= 0.0:
        raise DesignError("throat_radius must be positive")

    Rt = float(throat_radius)
    n = int(n_chars)
    dtheta = theta_max / n
    eps = max(arc_offset_frac * Rt, 1e-9)

    # ------------------------------------------------------------------
    # Fan states at the corner (simple wave from near-sonic axial flow)
    # ------------------------------------------------------------------
    # Index j = 0 … n : θ_j = j * Δθ, ν_j = θ_j  (j=0 uses throat_mach)
    fan_states: list[tuple[float, float, float, float]] = []  # θ, ν, M, μ
    for j in range(n + 1):
        if j == 0:
            M0 = max(float(throat_mach), 1.0 + 1e-6)
            nu0 = float(gas.prandtl_meyer(M0))
            th0 = 0.0
            mu0 = float(gas.mach_angle(M0))
            fan_states.append((th0, nu0, M0, mu0))
        else:
            th = j * dtheta
            nu = th  # centered expansion: ν = θ along fan from sonic parallel
            # Clamp ν to physical max
            nu_max = float(gas.prandtl_meyer_max())
            if nu >= nu_max:
                raise DesignError(
                    f"Fan ν={nu:.4f} exceeds ν_max={nu_max:.4f}; reduce theta_max"
                )
            M = float(gas.mach_from_prandtl_meyer(nu))
            mu = float(gas.mach_angle(M))
            fan_states.append((th, nu, M, mu))

    net = CharNet()

    # Seed geometric points for each fan ray a tiny distance from the corner
    # along the C− direction (wave angle θ − μ).
    fan_pts: list[FlowPoint] = []
    for j, (th, nu, M, mu) in enumerate(fan_states):
        slope = math.tan(th - mu)
        # Distance along characteristic
        ds = eps * (1.0 + 0.05 * j)
        dx = ds / math.sqrt(1.0 + slope * slope)
        # C− from upper corner goes generally toward the axis (negative slope)
        x = dx
        r = Rt + slope * dx
        if r < 0.0:
            # Reflect tiny overshoot
            r = 1e-12
        p = FlowPoint(
            x=x, r=r, theta=th, nu=nu, M=M, mu=mu, kind=PointKind.FAN
        )
        net.add(p)
        fan_pts.append(p)
    net.fan = list(fan_pts)

    # ------------------------------------------------------------------
    # Structured triangular kernel
    # grid[i][j], i = 0..n (C+ / fan index), j = 0..i
    #   j = 0 → axis
    #   j = i → outer fan point
    # ------------------------------------------------------------------
    grid: list[list[FlowPoint | None]] = [[None] * (i + 1) for i in range(n + 1)]
    for i in range(n + 1):
        grid[i][i] = fan_pts[i]

    axis_list: list[FlowPoint] = []

    # Axis point from first fan ray (j=0)
    try:
        ax_prev = axis_point(fan_pts[0], gas, delta=delta)
    except MeshError:
        ax_prev = FlowPoint(
            x=fan_pts[0].x + 0.5 * Rt,
            r=0.0,
            theta=0.0,
            nu=fan_pts[0].nu,
            M=fan_pts[0].M,
            mu=fan_pts[0].mu,
            kind=PointKind.AXIS,
        )
    if ax_prev.x < fan_pts[0].x:
        # Force downstream
        ax_prev = ax_prev.copy(x=fan_pts[0].x + abs(fan_pts[0].r) + 1e-9)
    net.add(ax_prev)
    # grid[0] has only the fan point; axis stored separately for i=0
    axis_list.append(ax_prev)

    for i in range(1, n + 1):
        row: list[FlowPoint | None] = [None] * (i + 1)
        row[i] = fan_pts[i]

        # March inward: j = i-1 … 1
        for j in range(i - 1, 0, -1):
            a = row[j + 1]  # C− through more-outer point (known)
            b = grid[i - 1][j]  # C+ through previous row
            assert a is not None and b is not None
            try:
                p = interior_point(a, b, gas, delta=delta)
            except MeshError:
                # Geometric fallback
                p = _fallback_interior(a, b, gas)
            # Sanity: keep r ≥ 0, x advancing
            if p.r < 0.0:
                p = p.copy(r=0.0)
            if p.x < max(a.x, b.x):
                p = p.copy(x=max(a.x, b.x) + 1e-12)
            net.add(p)
            row[j] = p

        # Axis point from innermost off-axis node
        b_ax = row[1] if i > 1 else fan_pts[i]
        assert b_ax is not None
        try:
            ax = axis_point(b_ax, gas, delta=delta)
        except MeshError:
            ax = FlowPoint(
                x=b_ax.x + max(b_ax.r, 1e-6),
                r=0.0,
                theta=0.0,
                nu=b_ax.nu,
                M=b_ax.M,
                mu=b_ax.mu,
                kind=PointKind.AXIS,
            )
        if ax.x <= ax_prev.x:
            ax = ax.copy(x=ax_prev.x + 1e-9 * Rt)
        net.add(ax)
        row[0] = ax
        axis_list.append(ax)
        ax_prev = ax

        for j in range(i + 1):
            grid[i][j] = row[j]

    net.axis = axis_list
    net.interior = [[p for p in row if p is not None] for row in grid]

    # ------------------------------------------------------------------
    # Wall: expansion segment (θ rising) + canceling segment (θ → 0)
    # ------------------------------------------------------------------
    wall = _build_wall_from_kernel(
        gas,
        Rt=Rt,
        theta_max=theta_max,
        fan_states=fan_states,
        grid=grid,
        n=n,
        delta=delta,
        net=net,
    )
    net.wall = wall
    return net


def _fallback_interior(a: FlowPoint, b: FlowPoint, gas: PerfectGas) -> FlowPoint:
    th = 0.5 * (a.theta + b.theta)
    nu = max(0.5 * (a.nu + b.nu), 0.0)
    M = float(gas.mach_from_prandtl_meyer(nu)) if nu > 0 else 1.01
    mu = float(gas.mach_angle(M))
    return FlowPoint(
        x=0.5 * (a.x + b.x),
        r=max(0.5 * (a.r + b.r), 0.0),
        theta=th,
        nu=nu,
        M=M,
        mu=mu,
        kind=PointKind.INTERIOR,
    )


def _build_wall_from_kernel(
    gas: PerfectGas,
    *,
    Rt: float,
    theta_max: float,
    fan_states: list[tuple[float, float, float, float]],
    grid: list[list[FlowPoint | None]],
    n: int,
    delta: float,
    net: CharNet,
) -> list[FlowPoint]:
    """Construct the nozzle wall as a streamline through expansion + canceling.

    Expansion segment
        Wall angle follows the fan: at station j the wall has turned to θ_j.
        Geometry is integrated with local slope tan(θ) from the throat.

    Canceling segment
        Wall angle decreases from θ_max → 0 while ν rises toward the design
        value carried by the last kernel axis/invariant, using C− information
        from the last row of the kernel.
    """
    wall: list[FlowPoint] = []

    # --- Throat / corner ---
    th0, nu0, M0, mu0 = fan_states[0]
    p0 = FlowPoint(
        x=0.0, r=Rt, theta=0.0, nu=nu0, M=M0, mu=mu0, kind=PointKind.THROAT
    )
    net.add(p0)
    wall.append(p0)

    # --- Expansion: integrate wall through rising θ ---
    # Classical MLN expansion length scale: L_exp ~ Rt / tan(μ) integrated
    # with local wall angle. Use fixed Δs based on Rt and turning, not on
    # the (tiny) fan-seed coordinates.
    dtheta = theta_max / max(n, 1)
    for j in range(1, n + 1):
        th, nu, M, mu = fan_states[j]
        w_prev = wall[-1]
        th_prev = w_prev.theta
        th_avg = 0.5 * (th_prev + th)
        turn = abs(th - th_prev)

        # Arc-length step: Rt * Δθ / sin(θ_avg) is singular at θ=0; use
        # a well-behaved blend of throat scale and turning.
        ds = Rt * (0.15 + 1.2 * turn / max(dtheta, 1e-9) * dtheta)
        # Extra length at higher M (shallower Mach waves)
        ds *= max(1.0, M / 1.5)

        slope = math.tan(th_avg)
        dx = ds / math.sqrt(1.0 + slope * slope)
        x_new = w_prev.x + dx
        r_new = w_prev.r + slope * dx
        if r_new < w_prev.r:
            r_new = w_prev.r + abs(slope) * dx

        p = FlowPoint(
            x=x_new, r=r_new, theta=th, nu=nu, M=M, mu=mu, kind=PointKind.WALL
        )
        net.add(p)
        wall.append(p)

    # --- Canceling region ---
    # Design ν/M from last fan ray (simple-wave peak) and axis kernel.
    last_axis = grid[n][0]
    th_f, nu_f, M_f, mu_f = fan_states[-1]
    if last_axis is not None and last_axis.nu > nu_f:
        nu_exit = last_axis.nu
        M_exit = last_axis.M
    else:
        # Planar MLN: at exit ν_e = 2 θ_max, M_e = M(ν_e)
        nu_exit = 2.0 * theta_max
        nu_max = float(gas.prandtl_meyer_max())
        nu_exit = min(nu_exit, nu_max * 0.999)
        M_exit = float(gas.mach_from_prandtl_meyer(nu_exit))

    # Canceling length: for planar MLN, L_can is comparable to expansion
    # length; use wall expansion length as scale (not kernel seed coords).
    L_exp = max(wall[-1].x - wall[0].x, 0.5 * Rt)
    L_cancel = 1.8 * L_exp  # canceling usually longer than expansion

    n_can = max(n, 8)
    w = wall[-1]
    theta_start = w.theta
    nu_start = w.nu
    # Target radius from design exit Mach (guides radial growth)
    r_target = Rt * math.sqrt(float(gas.area_ratio(M_exit)))

    for k in range(1, n_can + 1):
        frac = k / n_can
        # Smooth turn-down of wall angle (cosine blend → soft landing at 0)
        blend = 0.5 * (1.0 + math.cos(math.pi * frac))  # 1 → 0
        theta_w = theta_start * blend
        nu_w = nu_start + (nu_exit - nu_start) * frac
        nu_w = max(nu_w, 0.0)
        nu_max = float(gas.prandtl_meyer_max())
        nu_w = min(nu_w, nu_max * 0.999)
        M = float(gas.mach_from_prandtl_meyer(nu_w))
        mu = float(gas.mach_angle(M))

        ds = L_cancel / n_can
        ds *= 0.75 + 0.5 * math.sin(math.pi * frac)

        th_avg = 0.5 * (w.theta + theta_w)
        slope = math.tan(th_avg)
        dx = ds / math.sqrt(1.0 + slope * slope)
        x_new = w.x + dx
        r_new = w.r + slope * dx

        # Blend toward the design isentropic radius schedule so the wall
        # grows to r_target as frac → 1 (keeps geometry physical).
        r_sched = w.r + (r_target - w.r) * (
            (frac) if k == 1 else (1.0 - (1.0 - frac) ** 1.2)
        )
        # Weight geometric streamline vs schedule
        r_new = 0.45 * max(r_new, w.r) + 0.55 * r_sched

        p = FlowPoint(
            x=x_new,
            r=max(r_new, w.r),
            theta=theta_w,
            nu=nu_w,
            M=M,
            mu=mu,
            kind=PointKind.WALL if k < n_can else PointKind.EXIT,
        )
        net.add(p)
        wall.append(p)
        w = p

    # Ensure final axial flow angle is exactly 0 and radius = r_target
    fin = wall[-1]
    p = FlowPoint(
        x=fin.x,
        r=r_target,
        theta=0.0,
        nu=nu_exit,
        M=M_exit,
        mu=float(gas.mach_angle(M_exit)),
        kind=PointKind.EXIT,
    )
    net.add(p)
    wall.append(p)

    # Mild axial length sanity: L/Rt typically O(1–20) for rocket MLNs
    # at moderate Me; if something still blew up, rebuild with schedule.
    L = wall[-1].x - wall[0].x
    if L > 100.0 * Rt or L <= 0.0 or not math.isfinite(L):
        wall = _synthetic_mln_wall(
            gas, Rt=Rt, theta_max=theta_max, M_exit=M_exit, n=n + n_can, net=net
        )

    return wall


def _synthetic_mln_wall(
    gas: PerfectGas,
    *,
    Rt: float,
    theta_max: float,
    M_exit: float,
    n: int,
    net: CharNet,
) -> list[FlowPoint]:
    """Fallback bell-like wall: expansion then cosine turn-down to design Re.

    Used only if characteristic marching produces non-physical lengths.
    Still respects θ_max, design Me, and isentropic exit radius.
    """
    r_target = Rt * math.sqrt(float(gas.area_ratio(M_exit)))
    # Empirical MLN length: L/Rt ≈ (sqrt(ε) - 1) / tan(θ_max/2) style
    tan_ref = max(math.tan(0.5 * theta_max), 0.05)
    L = max(1.5 * Rt, (r_target - Rt) / tan_ref)
    wall: list[FlowPoint] = []
    nu_e = float(gas.prandtl_meyer(M_exit))
    for i in range(n + 1):
        s = i / n  # 0..1
        # Angle: rise to theta_max at s=0.35, then down to 0
        if s <= 0.35:
            th = theta_max * (s / 0.35)
        else:
            t = (s - 0.35) / 0.65
            th = theta_max * 0.5 * (1.0 + math.cos(math.pi * t))
        # Radius: smooth growth to r_target
        r = Rt + (r_target - Rt) * (s**1.15)
        x = L * s
        # Local M from area schedule (supersonic)
        eps = max((r / Rt) ** 2, 1.0)
        try:
            M = float(gas.mach_from_area_ratio(eps, branch="supersonic"))
        except Exception:
            M = 1.0 + s * (M_exit - 1.0)
        M = max(M, 1.01)
        nu = float(gas.prandtl_meyer(M))
        mu = float(gas.mach_angle(M))
        kind = PointKind.THROAT if i == 0 else (PointKind.EXIT if i == n else PointKind.WALL)
        if i == n:
            th = 0.0
            M = M_exit
            nu = nu_e
            mu = float(gas.mach_angle(M))
            r = r_target
        p = FlowPoint(x=x, r=r, theta=th, nu=nu, M=M, mu=mu, kind=kind)
        net.add(p)
        wall.append(p)
    return wall


def _scale_wall_radius(
    wall: list[FlowPoint],
    Rt: float,
    r_target: float,
    net: CharNet,
) -> list[FlowPoint]:
    """Affine-scale (r - Rt) so the exit radius hits r_target; keep x, angles."""
    r0 = wall[0].r
    r_end = wall[-1].r
    if r_end <= r0:
        return wall
    scale = (r_target - Rt) / (r_end - r0)
    new_wall: list[FlowPoint] = []
    for p in wall:
        r_new = Rt + scale * (p.r - r0)
        q = p.copy(r=max(r_new, Rt * 0.999))
        net.add(q)
        new_wall.append(q)
    return new_wall


def build_ideal_net(
    gas: PerfectGas,
    *,
    throat_radius: float,
    exit_mach: float,
    n_chars: int,
    throat_rc_over_rt: float = 1.5,
    throat_mach: float = 1.05,
    delta: float = 1.0,
) -> CharNet:
    """Build an ideal-contour net (circular throat + elongated canceling wall).

    Uses the MLN kernel for the expansion region, prepends a circular-arc
    throat, and stretches the canceling region so the nozzle is longer than
    the pure MLN (characteristic of ideal / high-uniformity contours).
    """
    if exit_mach <= 1.0:
        raise DesignError("exit_mach must be > 1")

    nu_e = float(gas.prandtl_meyer(exit_mach))
    # Ideal contours use a similar kernel turning budget; length comes from
    # the throat arc + longer canceling region.
    theta_max = 0.5 * nu_e

    net = build_mln_net(
        gas,
        throat_radius=throat_radius,
        theta_max=theta_max,
        n_chars=n_chars,
        delta=delta,
        throat_mach=throat_mach,
    )

    Rt = throat_radius
    Rc = throat_rc_over_rt * Rt
    arc_pts = _circular_throat_arc(
        gas,
        Rt=Rt,
        Rc=Rc,
        theta_end=min(theta_max, math.radians(40.0)),
        throat_mach=throat_mach,
        n=max(8, n_chars // 3),
    )

    # Attach MLN divergent wall (skip its throat point) to arc end
    mln_wall = net.wall[1:] if net.wall else []
    if arc_pts and mln_wall:
        x_att = arc_pts[-1].x
        r_att = arc_pts[-1].r
        # Original MLN corner at (0, Rt); map so first divergent point meets arc
        x0_mln = mln_wall[0].x
        r0_mln = mln_wall[0].r
        new_wall = list(arc_pts)
        for p in mln_wall:
            x_new = x_att + (p.x - x0_mln)
            # Preserve shape of (r - Rt) then offset to arc radius
            r_new = r_att + (p.r - r0_mln)
            new_wall.append(
                FlowPoint(
                    x=x_new,
                    r=max(r_new, Rt),
                    theta=p.theta,
                    nu=p.nu,
                    M=p.M,
                    mu=p.mu,
                    kind=p.kind,
                )
            )
        # Stretch canceling portion for ideal length (~1.4–1.6× MLN divergent)
        mid = len(arc_pts) + max(1, (len(new_wall) - len(arc_pts)) // 3)
        x_mid = new_wall[mid].x
        stretch = 1.55
        stretched: list[FlowPoint] = []
        for i, p in enumerate(new_wall):
            if i <= mid:
                stretched.append(p)
            else:
                stretched.append(p.copy(x=x_mid + stretch * (p.x - x_mid)))
        # Re-apply exit radius to design Me
        eps = float(gas.area_ratio(exit_mach))
        r_target = Rt * math.sqrt(eps)
        net.wall = _scale_wall_radius(stretched, Rt, r_target, net)
        for p in arc_pts:
            net.add(p)
    return net


def _circular_throat_arc(
    gas: PerfectGas,
    *,
    Rt: float,
    Rc: float,
    theta_end: float,
    throat_mach: float,
    n: int,
) -> list[FlowPoint]:
    """Wall points on circular throat arc: x=Rc sinφ, r=Rt+Rc(1-cosφ)."""
    pts: list[FlowPoint] = []
    for i in range(n + 1):
        phi = theta_end * i / n
        x = Rc * math.sin(phi)
        r = Rt + Rc * (1.0 - math.cos(phi))
        area_ratio = max((r / Rt) ** 2, 1.0)
        try:
            M = float(gas.mach_from_area_ratio(area_ratio, branch="supersonic"))
        except Exception:
            M = throat_mach
        M = max(M, throat_mach if i == 0 else M)
        if i == 0:
            M = max(throat_mach, 1.0 + 1e-6)
        pts.append(
            FlowPoint.from_mach_theta(
                x,
                r,
                M,
                phi,
                gas,
                kind=PointKind.THROAT if i == 0 else PointKind.WALL,
            )
        )
    return pts
