"""Sauer's transonic throat solution (small-perturbation, circular-arc throat).

Reference: Sauer (1947); reproduced in Mendelson, NASA TM X-53084 (1964),
Eqs. 2.5–2.6, and Zucrow & Hoffman Vol. 2. With δ = 0 (planar) or 1
(axisymmetric), Rt the throat radius and Rc the wall radius of curvature::

    α  = √((1 + δ) / ((γ + 1) Rc Rt))
    u  = α x + (γ + 1) α² r² / (2(1 + δ))
    v  = (γ + 1) α² x r / (1 + δ) + (γ + 1)² α³ r³ / (2(1 + δ)(3 + δ))

u, v are perturbations of V/a*; x is measured from the axial sonic point,
which lies ε = (γ + 1) α Rt² / (2(3 + δ)) downstream of the throat plane.
Sauer is accurate for Rc/Rt ≳ 2; Hall's higher-order series is not included.
"""

from __future__ import annotations

import math


def sauer_alpha(gamma: float, Rt: float, Rc: float, delta: float = 1.0) -> float:
    """Axial velocity gradient α at the sonic point [1/length]."""
    dl = 0.0 if abs(delta) < 0.5 else 1.0
    return math.sqrt((1.0 + dl) / ((gamma + 1.0) * Rc * Rt))


def sauer_offset(gamma: float, Rt: float, Rc: float, delta: float = 1.0) -> float:
    """Axial distance from the throat plane to the axial sonic point."""
    dl = 0.0 if abs(delta) < 0.5 else 1.0
    return (gamma + 1.0) * sauer_alpha(gamma, Rt, Rc, delta) * Rt * Rt / (2.0 * (3.0 + dl))


def sauer_uv(x_s: float, r: float, gamma: float, Rt: float, Rc: float, delta: float = 1.0) -> tuple[float, float]:
    """Perturbation velocities (u, v) at Sauer coordinate x_s (from sonic point) and radius r."""
    dl = 0.0 if abs(delta) < 0.5 else 1.0
    a = sauer_alpha(gamma, Rt, Rc, delta)
    gp = gamma + 1.0
    u = a * x_s + gp * a * a * r * r / (2.0 * (1.0 + dl))
    v = gp * a * a * x_s * r / (1.0 + dl) + gp * gp * a**3 * r**3 / (2.0 * (1.0 + dl) * (3.0 + dl))
    return u, v


def mach_from_mstar(ms: float, gamma: float) -> float:
    """M from the characteristic Mach number M* = V/a*."""
    g = gamma
    ms = min(ms, 0.999 * math.sqrt((g + 1.0) / (g - 1.0)))
    return math.sqrt(2.0 * ms * ms / ((g + 1.0) - (g - 1.0) * ms * ms))


def mstar_from_mach(M: float, gamma: float) -> float:
    g = gamma
    return math.sqrt(0.5 * (g + 1.0) * M * M / (1.0 + 0.5 * (g - 1.0) * M * M))


def sauer_line(
    gamma: float, Rt: float, Rc: float, *, axis_mach: float, n: int, delta: float = 1.0
) -> tuple[float, float, list[tuple[float, float, float]]]:
    """Vertical initial-data line where the axis Mach equals ``axis_mach``.

    Returns (x_line, φ_wall, [(r, M, θ), …]) with x measured from the throat
    plane, the wall point on the circular arc x = Rc sinφ, r = Rt + Rc(1 − cosφ).
    """
    a = sauer_alpha(gamma, Rt, Rc, delta)
    xs = (mstar_from_mach(axis_mach, gamma) - 1.0) / a
    x_line = xs + sauer_offset(gamma, Rt, Rc, delta)
    if x_line >= Rc:
        raise ValueError("axis_mach too high for this throat radius of curvature")
    phi = math.asin(x_line / Rc)
    r_w = Rt + Rc * (1.0 - math.cos(phi))
    pts = []
    for j in range(n + 1):
        r = r_w * j / n
        u, v = sauer_uv(xs, r, gamma, Rt, Rc, delta)
        th = phi if j == n else math.atan2(v, 1.0 + u)
        pts.append((r, mach_from_mstar(1.0 + u, gamma), th))
    return x_line, phi, pts
