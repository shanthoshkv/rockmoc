"""Unit processes: interior, axis, and wall points for axisymmetric MOC."""

from __future__ import annotations

import math

from rockmoc.exceptions import MeshError
from rockmoc.moc.equations import (
    average_state,
    char_slope_minus,
    char_slope_plus,
    source_minus,
    source_plus,
)
from rockmoc.moc.types import FlowPoint, PointKind
from rockmoc.thermo.gas_protocol import MOCGas


def _sync_mach(point_nu: float, gas: MOCGas) -> tuple[float, float, float]:
    """Return (M, μ, γ_local) using table-aware PM inverse when available."""
    M = float(gas.mach_from_prandtl_meyer(point_nu))
    mu = float(gas.mach_angle(M))
    g_loc = float(getattr(gas, "gamma_at_mach", lambda _m: gas.gamma)(M))
    return M, mu, g_loc


def interior_point(
    a: FlowPoint,
    b: FlowPoint,
    gas: MOCGas,
    *,
    delta: float = 1.0,
    max_iter: int = 12,
    tol: float = 1e-9,
) -> FlowPoint:
    """Compute interior point P from C− through A and C+ through B.

    A is the upstream point on the C− characteristic, B on the C+.
    """
    # Predictor: use properties at known points.
    lam_m = char_slope_minus(a.theta, a.mu)
    lam_p = char_slope_plus(b.theta, b.mu)
    denom = lam_m - lam_p
    if abs(denom) < 1e-16:
        raise MeshError("C+ and C− are parallel; cannot form interior point")

    x = (b.r - a.r - lam_p * b.x + lam_m * a.x) / denom
    r = a.r + lam_m * (x - a.x)
    if r < 0.0:
        r = 0.0

    sm = source_minus(a.theta, a.mu, a.r, delta=delta)
    sp = source_plus(b.theta, b.mu, b.r, delta=delta)
    km = (a.theta - a.nu) + sm * (x - a.x)  # θ − ν invariant + source
    kp = (b.theta + b.nu) + sp * (x - b.x)  # θ + ν invariant + source
    theta = 0.5 * (km + kp)
    nu = 0.5 * (kp - km)
    if nu < 0.0:
        nu = 0.0

    M, mu, g_loc = _sync_mach(nu, gas)
    p = FlowPoint(x=x, r=r, theta=theta, nu=nu, M=M, mu=mu, gamma_local=g_loc, kind=PointKind.INTERIOR)

    for _ in range(max_iter):
        th_am, mu_am, r_am, _nu_am = average_state(a, p)
        th_bp, mu_bp, r_bp, _nu_bp = average_state(b, p)
        # Use averaged radii carefully near axis
        r_am = max(r_am, 0.5 * (a.r + p.r))
        r_bp = max(r_bp, 0.5 * (b.r + p.r))

        lam_m = char_slope_minus(th_am, mu_am)
        lam_p = char_slope_plus(th_bp, mu_bp)
        denom = lam_m - lam_p
        if abs(denom) < 1e-16:
            raise MeshError("C+ and C− became parallel during corrector")

        x_new = (b.r - a.r - lam_p * b.x + lam_m * a.x) / denom
        r_new = a.r + lam_m * (x_new - a.x)
        if r_new < 0.0:
            r_new = 0.0

        # Under-relaxed average of predictor and corrector sources
        sm = source_minus(th_am, mu_am, max(r_am, 1e-14), delta=delta)
        sp = source_plus(th_bp, mu_bp, max(r_bp, 1e-14), delta=delta)
        km = (a.theta - a.nu) + sm * (x_new - a.x)
        kp = (b.theta + b.nu) + sp * (x_new - b.x)
        theta_new = 0.5 * (km + kp)
        nu_new = 0.5 * (kp - km)
        if nu_new < 0.0:
            nu_new = 0.0
        # Under-relaxation for stiff axisym sources
        omega = 0.65
        theta_new = (1.0 - omega) * p.theta + omega * theta_new
        nu_new = (1.0 - omega) * p.nu + omega * nu_new

        if (
            abs(x_new - p.x) < tol
            and abs(r_new - p.r) < tol
            and abs(theta_new - p.theta) < tol
            and abs(nu_new - p.nu) < tol
        ):
            M, mu, g_loc = _sync_mach(nu_new, gas)
            return FlowPoint(
                x=x_new,
                r=r_new,
                theta=theta_new,
                nu=nu_new,
                M=M,
                mu=mu,
                kind=PointKind.INTERIOR,
            )

        M, mu, g_loc = _sync_mach(nu_new, gas)
        p = FlowPoint(
            x=x_new,
            r=r_new,
            theta=theta_new,
            nu=nu_new,
            M=M,
            mu=mu,
            kind=PointKind.INTERIOR,
        )

    return p


def axis_point(
    b: FlowPoint,
    gas: MOCGas,
    *,
    delta: float = 1.0,
    max_iter: int = 12,
    tol: float = 1e-9,
) -> FlowPoint:
    """Axis point P (r=0, θ=0) reached by C+ through off-axis point B.

    Symmetry: the C− through P is the mirror of C+, so only one independent
    compatibility equation is required. On the axis the geometric source is
    treated in the limit (vanishes for θ→0).
    """
    # C+ from B: r_P = 0 = r_B + tan(θ_B+μ_B)*(x_P - x_B)
    lam_p = char_slope_plus(b.theta, b.mu)
    if abs(lam_p) < 1e-16:
        raise MeshError("Cannot reach axis: C+ is horizontal")
    x = b.x - b.r / lam_p
    r = 0.0
    theta = 0.0

    # On axis, S+ → 0 when θ→0. Use planar-like invariant along C+ then correct.
    # (θ + ν)_P = (θ + ν)_B + S+_avg*(x_P - x_B), θ_P=0
    sp = source_plus(b.theta, b.mu, b.r, delta=delta)
    nu = (b.theta + b.nu) + sp * (x - b.x)  # since θ_P=0 → ν_P = (θ+ν)_P
    if nu < 0.0:
        nu = 0.0
    M, mu, g_loc = _sync_mach(nu, gas)
    p = FlowPoint(x=x, r=r, theta=theta, nu=nu, M=M, mu=mu, gamma_local=g_loc, kind=PointKind.AXIS)

    for _ in range(max_iter):
        th_bp, mu_bp, r_bp, _ = average_state(b, p)
        r_bp = max(0.5 * (b.r + p.r), 1e-14)
        lam_p = char_slope_plus(th_bp, mu_bp)
        if abs(lam_p) < 1e-16:
            raise MeshError("C+ horizontal during axis corrector")
        x_new = b.x - b.r / lam_p
        sp = source_plus(th_bp, mu_bp, r_bp, delta=delta)
        # Limit source as we approach axis: blend to 0
        sp *= b.r / (b.r + p.r + 1e-14)
        nu_new = (b.theta + b.nu) + sp * (x_new - b.x)
        # θ_P = 0
        if nu_new < 0.0:
            nu_new = 0.0

        if abs(x_new - p.x) < tol and abs(nu_new - p.nu) < tol:
            M, mu, g_loc = _sync_mach(nu_new, gas)
            return FlowPoint(
                x=x_new, r=0.0, theta=0.0, nu=nu_new, M=M, mu=mu, gamma_local=g_loc, kind=PointKind.AXIS
            )

        M, mu, g_loc = _sync_mach(nu_new, gas)
        p = FlowPoint(
            x=x_new, r=0.0, theta=0.0, nu=nu_new, M=M, mu=mu, gamma_local=g_loc, kind=PointKind.AXIS
        )

    return p


def wall_point(
    a: FlowPoint,
    wall_theta: float,
    gas: MOCGas,
    *,
    delta: float = 1.0,
    max_iter: int = 40,
    tol: float = 1e-10,
) -> FlowPoint:
    """Wall point P on a prescribed wall inclination ``wall_theta``.

    C− through known interior/fan point A intersects the wall streamline
    direction ``dr/dx = tan(wall_theta)`` (locally constant wall angle).
    """
    # Wall line through previous wall point is handled by caller via a being
    # the C− point and wall passing through a reference; here we assume the
    # wall slope is wall_theta and the wall passes through... actually for
    # successive wall points, A is interior (C−) and we also need previous
    # wall point W. This function uses a pure C− + known wall angle form
    # with position closed by intersection with line from a reference.
    raise MeshError("Use wall_point_from_previous for sequential wall marching")


def wall_point_from_previous(
    a: FlowPoint,
    w_prev: FlowPoint,
    gas: MOCGas,
    *,
    wall_theta: float | None = None,
    delta: float = 1.0,
    max_iter: int = 40,
    tol: float = 1e-10,
) -> FlowPoint:
    """Next wall point from C− through A and wall streamline through W_prev.

    If ``wall_theta`` is None, the wall angle is unknown and determined by
    the canceling condition for a simple-wave residual: for pure canceling
    of the remaining C+ waves, θ_wall = θ along the last interior relation.

    Standard **canceling wall** for MLN: the wall is a streamline, and the
    flow angle at the wall equals the wall geometric angle. The C− from the
    interior determines (θ, ν) with the streamline constraint
    ``dr/dx = tan(θ)``.
    """
    # Unknown wall angle: use θ at the wall equal to flow angle (streamline).
    # Unit process: C− from A + streamline from W_prev with slope tan(θ_P).
    theta_guess = a.theta if wall_theta is None else wall_theta
    nu_guess = a.nu
    M, mu, g_loc = _sync_mach(max(nu_guess, 0.0), gas)
    p = FlowPoint(
        x=w_prev.x + 1e-6,
        r=w_prev.r,
        theta=theta_guess,
        nu=nu_guess,
        M=M,
        mu=mu,
        kind=PointKind.WALL,
    )

    for _ in range(max_iter):
        th_am, mu_am, r_am, _ = average_state(a, p)
        r_am = max(r_am, 1e-14)

        lam_m = char_slope_minus(th_am, mu_am)
        # Streamline slope uses wall flow angle
        if wall_theta is None:
            lam_w = math.tan(p.theta)
        else:
            lam_w = math.tan(wall_theta)

        denom = lam_m - lam_w
        if abs(denom) < 1e-16:
            # Nearly parallel: step in x slightly using C− only
            dx = max(1e-6 * (abs(a.x) + 1.0), 1e-8)
            x_new = a.x + dx
            r_new = a.r + lam_m * dx
        else:
            x_new = (w_prev.r - a.r - lam_w * w_prev.x + lam_m * a.x) / denom
            r_new = a.r + lam_m * (x_new - a.x)

        if r_new < 0.0:
            r_new = 0.0

        sm = source_minus(th_am, mu_am, r_am, delta=delta)
        # (θ − ν)_P = (θ − ν)_A + S−*(x_P - x_A)
        km = (a.theta - a.nu) + sm * (x_new - a.x)

        if wall_theta is None:
            # Streamline: flow angle is free but wall is streamline.
            # Second relation from remaining simple-wave intent is not
            # unique without C+; for canceling wall of MLN the standard
            # approach uses the fact that along the wall after the last
            # expansion, C+ invariants determine ν from θ.
            # Use: on canceling characteristic network, wall point also
            # satisfies C+ from a neighbor when available. Here A alone:
            # assume dθ = dν along wall for last simple wave cancel
            # → θ + ν = const from corner last ray... For sequential
            # interior-fed wall: θ_P is unknown, need another equation.
            #
            # Correct unit process for wall as streamline with C− only
            # requires known wall shape. For design (unknown wall), the
            # wall is determined by making the last family cancel:
            # θ_P = ν_e - ν_P  (planar) generalized with axisym sources.
            #
            # Engineering MLN wall construction: each wall point is the
            # streamline through previous wall with C− from interior
            # point on the last fan/characteristic of that "layer", and
            # θ_wall = θ_flow from solving C− with the design constraint
            # that residual turning equals remaining ν to reach axial
            # exit. For the classic construction, interior net is built
            # first; wall is the streamline of the outermost points.
            #
            # Practical approach used here: wall points track the
            # outermost characteristic points' streamline — caller should
            # pass wall_theta = theta of the design wall at this station
            # OR we set θ_P from C− by also using invariant from fan.
            #
            # For design wall (unknown shape), use:
            #   θ_P + ν_P = θ_last + ν_last  (along C+ from outer fan)
            # when a C+ partner is not supplied, fall back to θ_P = a.theta
            # for a first-order wall march of the outer streamline.
            theta_new = a.theta
            nu_new = theta_new - km  # from km = θ − ν
        else:
            theta_new = wall_theta
            nu_new = theta_new - km

        if nu_new < 0.0:
            nu_new = 0.0

        if (
            abs(x_new - p.x) < tol
            and abs(r_new - p.r) < tol
            and abs(theta_new - p.theta) < tol
            and abs(nu_new - p.nu) < tol
        ):
            M, mu, g_loc = _sync_mach(nu_new, gas)
            return FlowPoint(
                x=x_new,
                r=r_new,
                theta=theta_new,
                nu=nu_new,
                M=M,
                mu=mu,
                kind=PointKind.WALL,
            )

        M, mu, g_loc = _sync_mach(nu_new, gas)
        p = FlowPoint(
            x=x_new,
            r=r_new,
            theta=theta_new,
            nu=nu_new,
            M=M,
            mu=mu,
            kind=PointKind.WALL,
        )

    return p


def wall_point_canceling(
    a: FlowPoint,
    w_prev: FlowPoint,
    b_plus: FlowPoint,
    gas: MOCGas,
    *,
    delta: float = 1.0,
    max_iter: int = 40,
    tol: float = 1e-10,
) -> FlowPoint:
    """Wall point using C− from A, C+ from B_plus, constrained to streamline.

    The wall location is the intersection of the streamline through
    ``w_prev`` (slope tan θ_P) with the characteristic solution from A and
    B. This is the standard design unit process for an unknown wall.
    """
    # Initial guess from interior of A and B
    p = interior_point(a, b_plus, gas, delta=delta, max_iter=max_iter, tol=tol)
    p.kind = PointKind.WALL

    for _ in range(max_iter):
        # Streamline from w_prev with slope tan(θ_P)
        lam_w = math.tan(p.theta)
        th_am, mu_am, r_am, _ = average_state(a, p)
        th_bp, mu_bp, r_bp, _ = average_state(b_plus, p)
        r_am = max(r_am, 1e-14)
        r_bp = max(r_bp, 1e-14)

        lam_m = char_slope_minus(th_am, mu_am)
        lam_p = char_slope_plus(th_bp, mu_bp)

        # Intersect C− from A with streamline from w_prev
        denom = lam_m - lam_w
        if abs(denom) < 1e-16:
            x_new = p.x
            r_new = p.r
        else:
            x_new = (w_prev.r - a.r - lam_w * w_prev.x + lam_m * a.x) / denom
            r_new = a.r + lam_m * (x_new - a.x)
        if r_new < 0.0:
            r_new = abs(r_new)

        sm = source_minus(th_am, mu_am, r_am, delta=delta)
        sp = source_plus(th_bp, mu_bp, r_bp, delta=delta)
        km = (a.theta - a.nu) + sm * (x_new - a.x)
        kp = (b_plus.theta + b_plus.nu) + sp * (x_new - b_plus.x)
        theta_new = 0.5 * (km + kp)
        nu_new = 0.5 * (kp - km)
        if nu_new < 0.0:
            nu_new = 0.0

        # Recompute streamline slope consistency: wall angle = flow angle
        lam_w2 = math.tan(theta_new)
        denom2 = lam_m - lam_w2
        if abs(denom2) > 1e-16:
            x_new = (w_prev.r - a.r - lam_w2 * w_prev.x + lam_m * a.x) / denom2
            r_new = a.r + lam_m * (x_new - a.x)

        if (
            abs(x_new - p.x) < tol
            and abs(r_new - p.r) < tol
            and abs(theta_new - p.theta) < tol
            and abs(nu_new - p.nu) < tol
        ):
            M, mu, g_loc = _sync_mach(nu_new, gas)
            return FlowPoint(
                x=x_new,
                r=max(r_new, 0.0),
                theta=theta_new,
                nu=nu_new,
                M=M,
                mu=mu,
                kind=PointKind.WALL,
            )

        M, mu, g_loc = _sync_mach(nu_new, gas)
        p = FlowPoint(
            x=x_new,
            r=max(r_new, 0.0),
            theta=theta_new,
            nu=nu_new,
            M=M,
            mu=mu,
            kind=PointKind.WALL,
        )
        # silence unused lam_p (kept for future C+ intersection variant)
        _ = lam_p

    return p
