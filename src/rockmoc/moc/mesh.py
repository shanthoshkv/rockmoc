"""Unit processes: interior, axis, and wall points for axisymmetric MOC."""

from __future__ import annotations

from rockmoc import _kernels as _core
from rockmoc.exceptions import MeshError, ThermoError
from rockmoc.moc.equations import (
    average_state,
    char_slope_minus,
    char_slope_plus,
    source_minus,
    source_plus,
)
from rockmoc.moc.types import FlowPoint, PointKind
from rockmoc.thermo.gas_protocol import MOCGas


def _fast_point(x: float, r: float, th: float, nu: float, g: float, kind: PointKind) -> FlowPoint:
    if nu > _core.pm_max(g) + 1e-12:
        raise ThermoError(f"ν exceeds maximum for γ={g}")
    M = _core.pm_inv(nu, g)
    return FlowPoint(x=x, r=r, theta=th, nu=nu, M=M, mu=_core.mach_angle(M), kind=kind, gamma_local=g)


def _sync_mach(point_nu: float, gas: MOCGas) -> tuple[float, float, float]:
    """Return (M, μ, γ_local) using table-aware PM inverse when available."""
    M = float(gas.mach_from_prandtl_meyer(point_nu))
    mu = float(gas.mach_angle(M))
    g_loc = float(gas.gamma_at_mach(M))
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
    g = _core.const_gamma(gas)
    if g is not None:
        x, r, th, nu, ok = _core.interior(
            a.x, a.r, a.theta, a.nu, a.mu, b.x, b.r, b.theta, b.nu, b.mu,
            g, float(delta), max_iter, tol,
        )
        if not ok:
            raise MeshError("C+ and C− are parallel; cannot form interior point")
        return _fast_point(x, r, th, nu, g, PointKind.INTERIOR)
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
    km = (a.theta + a.nu) + sm * (x - a.x)  # θ + ν along C− + source
    kp = (b.theta - b.nu) + sp * (x - b.x)  # θ − ν along C+ + source
    theta = 0.5 * (km + kp)
    nu = 0.5 * (km - kp)
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
        km = (a.theta + a.nu) + sm * (x_new - a.x)
        kp = (b.theta - b.nu) + sp * (x_new - b.x)
        theta_new = 0.5 * (km + kp)
        nu_new = 0.5 * (km - kp)
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
    """Axis point P (r=0, θ=0) reached by the C− through off-axis point B.

    Symmetry: the C+ through P is the mirror of the C−, so only one independent
    compatibility equation is required: (θ+ν)_P = (θ+ν)_B + S−·Δx.
    """
    g = _core.const_gamma(gas)
    if g is not None:
        x, nu, ok = _core.axis(b.x, b.r, b.theta, b.nu, b.mu, g, float(delta), max_iter, tol)
        if not ok:
            raise MeshError("Cannot reach axis: C− is horizontal")
        return _fast_point(x, 0.0, 0.0, nu, g, PointKind.AXIS)
    # C− from B: r_P = 0 = r_B + tan(θ_B−μ_B)*(x_P − x_B)
    lam_p = char_slope_minus(b.theta, b.mu)
    if abs(lam_p) < 1e-16:
        raise MeshError("Cannot reach axis: C− is horizontal")
    x = b.x - b.r / lam_p
    r = 0.0
    theta = 0.0

    # (θ + ν)_P = (θ + ν)_B + S−_avg*(x_P - x_B), θ_P=0
    sp = source_minus(b.theta, b.mu, b.r, delta=delta)
    nu = (b.theta + b.nu) + sp * (x - b.x)  # since θ_P=0 → ν_P = (θ+ν)_P
    if nu < 0.0:
        nu = 0.0
    M, mu, g_loc = _sync_mach(nu, gas)
    p = FlowPoint(x=x, r=r, theta=theta, nu=nu, M=M, mu=mu, gamma_local=g_loc, kind=PointKind.AXIS)

    for _ in range(max_iter):
        th_bp, mu_bp, r_bp, _ = average_state(b, p)
        r_bp = max(0.5 * (b.r + p.r), 1e-14)
        lam_p = char_slope_minus(th_bp, mu_bp)
        if abs(lam_p) < 1e-16:
            raise MeshError("C− horizontal during axis corrector")
        x_new = b.x - b.r / lam_p
        sp = source_minus(th_bp, mu_bp, r_bp, delta=delta)
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
