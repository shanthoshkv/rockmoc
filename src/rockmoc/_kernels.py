"""Float-only MOC kernels for constant-γ gases (optionally Numba-compiled).

These mirror :mod:`rockmoc.moc.equations` / :mod:`rockmoc.moc.mesh` exactly but
work on plain floats, avoiding NumPy scalar overhead (~50× faster in CPython).
With ``numba`` installed (``pip install rockmoc[fast]``) they are JIT-compiled.
"""

from __future__ import annotations

import math

try:  # pragma: no cover - exercised only when numba is installed
    from numba import njit as _njit

    def jit(fn):  # type: ignore[no-untyped-def]
        return _njit(cache=True, fastmath=False)(fn)

    HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    def jit(fn):  # type: ignore[no-untyped-def]
        return fn

    HAVE_NUMBA = False


@jit
def pm(M: float, g: float) -> float:
    """Prandtl–Meyer ν(M) [rad]."""
    gp = math.sqrt((g + 1.0) / (g - 1.0))
    rt = math.sqrt(max(M * M - 1.0, 0.0))
    return gp * math.atan(rt / gp) - math.atan(rt)


@jit
def pm_max(g: float) -> float:
    return 0.5 * math.pi * (math.sqrt((g + 1.0) / (g - 1.0)) - 1.0)


@jit
def pm_inv(nu: float, g: float) -> float:
    """Inverse Prandtl–Meyer: safeguarded Newton with analytic dν/dM."""
    if nu <= 1e-15:
        return 1.0
    lo, hi = 1.0, 50.0
    while pm(hi, g) - nu < 0.0 and hi < 1e4:
        hi *= 1.5
    M = min(2.0, 0.5 * (lo + hi))
    for _ in range(100):
        fm = pm(M, g) - nu
        if fm > 0.0:
            hi = M
        else:
            lo = M
        if abs(fm) < 1e-15 or hi - lo < 1e-15 * hi:
            break
        d = math.sqrt(M * M - 1.0) / (M * (1.0 + 0.5 * (g - 1.0) * M * M))
        M_new = M - fm / d if d > 0.0 else 0.5 * (lo + hi)
        M = M_new if lo < M_new < hi else 0.5 * (lo + hi)
    return M


@jit
def area_ratio(M: float, g: float) -> float:
    expo = 0.5 * (g + 1.0) / (g - 1.0)
    return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** expo


@jit
def mach_angle(M: float) -> float:
    return math.asin(min(1.0, 1.0 / max(M, 1.0)))


@jit
def _source_core(theta: float, mu: float, r: float, phi: float, delta: float) -> float:
    if abs(delta) < 1e-14:
        return 0.0
    r_eff = r if abs(r) > 1e-10 else 1e-10
    s_over_r = theta / r_eff if abs(theta) < 1e-8 else math.sin(theta) / r_eff
    c = math.cos(phi)
    if abs(c) < 1e-12:
        return 0.0
    return -delta * s_over_r * math.sin(mu) / c


@jit
def source_plus(theta: float, mu: float, r: float, delta: float) -> float:
    """d(θ−ν)/dx along C+."""
    return _source_core(theta, mu, r, theta + mu, delta)


@jit
def source_minus(theta: float, mu: float, r: float, delta: float) -> float:
    """d(θ+ν)/dx along C−."""
    return -_source_core(theta, mu, r, theta - mu, delta)


@jit
def _nu_state(nu: float, g: float) -> tuple[float, float]:
    M = pm_inv(nu, g)
    return M, mach_angle(M)


@jit
def interior(
    xa: float, ra: float, tha: float, nua: float, mua: float,
    xb: float, rb: float, thb: float, nub: float, mub: float,
    g: float, delta: float, max_iter: int, tol: float,
) -> tuple[float, float, float, float, bool]:
    """Interior point from A (on C−) and B (on C+). Returns (x, r, θ, ν, ok)."""
    lam_m = math.tan(tha - mua)
    lam_p = math.tan(thb + mub)
    den = lam_m - lam_p
    if abs(den) < 1e-16:
        return 0.0, 0.0, 0.0, 0.0, False
    x = (rb - ra - lam_p * xb + lam_m * xa) / den
    r = max(ra + lam_m * (x - xa), 0.0)
    km = (tha + nua) + source_minus(tha, mua, ra, delta) * (x - xa)
    kp = (thb - nub) + source_plus(thb, mub, rb, delta) * (x - xb)
    th = 0.5 * (km + kp)
    nu = max(0.5 * (km - kp), 0.0)
    _, mu = _nu_state(nu, g)
    for _ in range(max_iter):
        th_am, mu_am, r_am = 0.5 * (tha + th), 0.5 * (mua + mu), 0.5 * (ra + r)
        th_bp, mu_bp, r_bp = 0.5 * (thb + th), 0.5 * (mub + mu), 0.5 * (rb + r)
        lam_m = math.tan(th_am - mu_am)
        lam_p = math.tan(th_bp + mu_bp)
        den = lam_m - lam_p
        if abs(den) < 1e-16:
            return 0.0, 0.0, 0.0, 0.0, False
        x_new = (rb - ra - lam_p * xb + lam_m * xa) / den
        r_new = max(ra + lam_m * (x_new - xa), 0.0)
        km = (tha + nua) + source_minus(th_am, mu_am, max(r_am, 1e-14), delta) * (x_new - xa)
        kp = (thb - nub) + source_plus(th_bp, mu_bp, max(r_bp, 1e-14), delta) * (x_new - xb)
        th_new = 0.65 * (0.5 * (km + kp)) + 0.35 * th
        nu_new = 0.65 * max(0.5 * (km - kp), 0.0) + 0.35 * nu
        done = (abs(x_new - x) < tol and abs(r_new - r) < tol
                and abs(th_new - th) < tol and abs(nu_new - nu) < tol)
        x, r, th, nu = x_new, r_new, th_new, nu_new
        _, mu = _nu_state(nu, g)
        if done:
            break
    return x, r, th, nu, True


@jit
def axis(
    xb: float, rb: float, thb: float, nub: float, mub: float,
    g: float, delta: float, max_iter: int, tol: float,
) -> tuple[float, float, bool]:
    """Axis point reached by the C− through B. Returns (x, ν, ok)."""
    lam = math.tan(thb - mub)
    if abs(lam) < 1e-16:
        return 0.0, 0.0, False
    x = xb - rb / lam
    nu = max((thb + nub) + source_minus(thb, mub, rb, delta) * (x - xb), 0.0)
    _, mu = _nu_state(nu, g)
    for _ in range(max_iter):
        th_bp, mu_bp = 0.5 * thb, 0.5 * (mub + mu)
        r_bp = max(0.5 * rb, 1e-14)
        lam = math.tan(th_bp - mu_bp)
        if abs(lam) < 1e-16:
            return 0.0, 0.0, False
        x_new = xb - rb / lam
        nu_new = max((thb + nub) + source_minus(th_bp, mu_bp, r_bp, delta) * (x_new - xb), 0.0)
        done = abs(x_new - x) < tol and abs(nu_new - nu) < tol
        x, nu = x_new, nu_new
        _, mu = _nu_state(nu, g)
        if done:
            break
    return x, nu, True


def const_gamma(gas: object) -> float | None:
    """γ if ``gas`` is a constant-γ perfect gas (fast path eligible), else None."""
    from rockmoc.thermo.perfect_gas import PerfectGas

    return gas.gamma if type(gas) is PerfectGas else None  # type: ignore[attr-defined]
