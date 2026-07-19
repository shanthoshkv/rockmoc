"""Robust scalar root finders used by thermo and MOC modules."""

from __future__ import annotations

from collections.abc import Callable

from rockmoc.exceptions import ConvergenceError


def bracketed_newton(
    f: Callable[[float], float],
    df: Callable[[float], float] | None,
    x0: float,
    *,
    x_min: float,
    x_max: float,
    tol: float = 1e-12,
    max_iter: int = 80,
    f_tol: float = 1e-14,
) -> float:
    """Solve f(x)=0 with Newton steps clamped to [x_min, x_max].

    Falls back to bisection-style steps when the Newton update leaves the
    bracket or the derivative is unusable.
    """
    if not (x_min < x_max):
        raise ConvergenceError(f"Invalid bracket [{x_min}, {x_max}]")

    x = min(max(x0, x_min), x_max)
    lo, hi = x_min, x_max
    f_lo = f(lo)
    f_hi = f(hi)
    if f_lo == 0.0:
        return lo
    if f_hi == 0.0:
        return hi
    if f_lo * f_hi > 0.0:
        raise ConvergenceError(
            f"Root not bracketed: f({lo})={f_lo}, f({hi})={f_hi}"
        )

    for _ in range(max_iter):
        fx = f(x)
        if abs(fx) < f_tol:
            return x

        # Maintain a bracketing interval when possible.
        if f_lo * fx <= 0.0:
            hi, f_hi = x, fx
        else:
            lo, f_lo = x, fx

        step_ok = False
        if df is not None:
            d = df(x)
            if abs(d) > 1e-16:
                x_new = x - fx / d
                if lo <= x_new <= hi:
                    x = x_new
                    step_ok = True

        if not step_ok:
            x = 0.5 * (lo + hi)

        if abs(hi - lo) < tol:
            return 0.5 * (lo + hi)

    raise ConvergenceError(
        f"Root finder failed after {max_iter} iterations; last x={x}, |f|={abs(f(x))}"
    )
