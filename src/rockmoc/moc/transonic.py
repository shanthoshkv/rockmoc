"""Sauer's transonic throat solution (small-perturbation, circular-arc throat).

Reference: Sauer (1947); reproduced in Mendelson, NASA TM X-53084 (1964),
Eqs. 2.5–2.6, and Zucrow & Hoffman Vol. 2. With δ = 0 (planar) or 1
(axisymmetric), Rt the throat radius and Rc the wall radius of curvature::

    α  = √((1 + δ) / ((γ + 1) Rc Rt))
    u  = α x + (γ + 1) α² r² / (2(1 + δ))
    v  = (γ + 1) α² x r / (1 + δ) + (γ + 1)² α³ r³ / (2(1 + δ)(3 + δ))

u, v are perturbations of V/a*; x is measured from the axial sonic point,
which lies ε = (γ + 1) α Rt² / (2(3 + δ)) downstream of the throat plane.
Sauer is accurate for Rc/Rt ≳ 2.

Hall (1962), axisymmetric, second order in ε = Rt/Rc (lengths / Rt, x from
the throat plane, z = x √(2R/(γ+1)), R = Rc/Rt)::

    u = 1 + ε u₁ + ε² u₂ ,   v = √((γ+1)ε/2) (ε v₁ + ε² v₂)
    u₁ = r²/2 − 1/4 + z
    v₁ = r³/4 − r/4 + r z
    u₂ = (2γ+9)/24 r⁴ − (4γ+15)/24 r² + (10γ+57)/288 + z (r² − 5/8) − (2γ−3)/6 z²
    v₂ = (γ+3)/9 r⁵ − (20γ+63)/96 r³ + (28γ+93)/288 r
         + z ((2γ+9)/6 r³ − (4γ+15)/12 r) + r z²

The first-order terms are Sauer's. Residuals of the exact potential equation
fall as ε^3.4 against ε^2.5 for Sauer (``tests/test_accuracy_v08.py``), but
the series diverges for tight throats: throat-plane axis Mach 0.93 at
Rc/Rt = 1.5, 1.0 at 1.0 and 1.30 at 0.625 (measured 0.8, Cuffel et al. 1969).
Hall is therefore used for Rc/Rt ≥ :data:`HALL_MIN_RC_RT` (it matches the
Rc/Rt = 2 wall pressure of NASA TN D-5377 to 0.3 % vs 5 % for Sauer).

Kliegel & Levine (AIAA J. 7(7), 1969) re-expand Hall's third-order series in
ε = 1/(R + 1) (substituting 1/R = ε/(1 − ε)), which stays bounded for tight
throats::

    u = 1 + ε u₁ + ε²(u₁ + u₂) + ε³(u₁ + 2u₂ + u₃)
    v = √((γ+1)ε/2) [ε v₁ + ε²(3/2 v₁ + v₂) + ε³(15/8 v₁ + 5/2 v₂ + v₃)]

with Hall's z = x √(2R/(γ+1)) (the only scaling for which ∂u/∂r = ∂v/∂x holds
order by order). u₃, v₃ are Hall's third-order terms as transcribed by Kbab
(Blida, 2012) and checked for irrotationality. At Rc/Rt = 0.625 it gives
throat-plane axis/wall Mach 0.79/1.31 and Cd 0.982 against 0.8/1.4 and 0.985
measured (Cuffel et al. 1969); Levine & Coats (NASA CR, 1970) report agreement
down to Rc/Rt ≈ 0.4. Used for :data:`KL_MIN_RC_RT` ≤ Rc/Rt < 1.5. Planar Hall
and KL are not included (Sauer used).
"""

from __future__ import annotations

import math

HALL_MIN_RC_RT = 1.5  # below this the Hall series diverges (see module docstring)
KL_MIN_RC_RT = 0.4  # Levine & Coats 1970: KL agrees with data down to ≈ 0.4


def _hall_terms(r: float, z: float, g: float) -> tuple[float, float, float, float, float, float]:
    """Hall's u₁…u₃, v₁…v₃ (lengths / Rt)."""
    r2, G = r * r, g * g
    u1 = 0.5 * r2 - 0.25 + z
    v1 = 0.25 * r2 * r - 0.25 * r + r * z
    u2 = ((2 * g + 9) / 24 * r2 * r2 - (4 * g + 15) / 24 * r2 + (10 * g + 57) / 288
          + z * (r2 - 0.625) - (2 * g - 3) / 6 * z * z)
    v2 = ((g + 3) / 9 * r**5 - (20 * g + 63) / 96 * r**3 + (28 * g + 93) / 288 * r
          + z * ((2 * g + 9) / 6 * r**3 - (4 * g + 15) / 12 * r) + r * z * z)
    u3 = ((556 * G + 1737 * g + 3069) / 10368 * r**6 - (388 * G + 1161 * g + 1881) / 2304 * r**4
          + (304 * G + 831 * g + 1242) / 1728 * r2 - (2708 * G + 7839 * g + 14211) / 82944
          + z * ((52 * G + 51 * g + 327) / 384 * r**4 - (52 * G + 75 * g + 279) / 192 * r2
                 + (92 * G + 180 * g + 639) / 1152)
          + z * z * (-(7 * g - 3) / 8 * r2 + (13 * g - 27) / 144) + (4 * G - 57 * g + 27) / 144 * z**3)
    v3 = ((6836 * G + 23031 * g + 30627) / 82944 * r**7 - (3380 * G + 11391 * g + 15291) / 13824 * r**5
          + (3424 * G + 11271 * g + 15228) / 13824 * r**3 - (7100 * G + 22311 * g + 30249) / 82944 * r
          + z * ((556 * G + 1737 * g + 3069) / 1728 * r**5 - (388 * G + 1161 * g + 1881) / 576 * r**3
                 + (304 * G + 831 * g + 1242) / 864 * r)
          + z * z * ((52 * G + 51 * g + 327) / 192 * r**3 - (52 * G + 75 * g + 279) / 192 * r)
          - (7 * g - 3) / 12 * r * z**3)
    return u1, u2, u3, v1, v2, v3


def kl_uv(x: float, r: float, gamma: float, Rt: float, Rc: float) -> tuple[float, float]:
    """Kliegel–Levine (u, v) = V/a* components; x from the throat plane."""
    R = Rc / Rt
    e = 1.0 / (R + 1.0)
    z = (x / Rt) * math.sqrt(2.0 * R / (gamma + 1.0))
    u1, u2, u3, v1, v2, v3 = _hall_terms(r / Rt, z, gamma)
    u = 1.0 + e * u1 + e * e * (u1 + u2) + e**3 * (u1 + 2.0 * u2 + u3)
    v = math.sqrt(0.5 * (gamma + 1.0) * e) * (e * v1 + e * e * (1.5 * v1 + v2)
                                              + e**3 * (1.875 * v1 + 2.5 * v2 + v3))
    return u, v


def discharge_coefficient(gamma: float, rc_over_rt: float) -> float:
    """Inviscid axisymmetric throat discharge coefficient, Kliegel–Levine series.

    1 − Cd = (γ+1)/(1+R)² [1/96 − (8γ−27)/(2304(1+R)) + (754γ²−757γ+3633)/(276480(1+R)²)]
    (converges for all R; agrees with the integrated KL velocity field to 5e-4).
    """
    g, s = gamma, 1.0 + rc_over_rt
    return 1.0 - (g + 1.0) / s**2 * (1.0 / 96.0 - (8.0 * g - 27.0) / (2304.0 * s)
                                     + (754.0 * g * g - 757.0 * g + 3633.0) / (276480.0 * s * s))


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


def hall_uv(x: float, r: float, gamma: float, Rt: float, Rc: float) -> tuple[float, float]:
    """Hall 2nd-order axisymmetric (u, v) = V/a* components; x from the throat plane."""
    R = Rc / Rt
    e = 1.0 / R
    z = (x / Rt) * math.sqrt(2.0 * R / (gamma + 1.0))
    u1, u2, _, v1, v2, _ = _hall_terms(r / Rt, z, gamma)
    u = 1.0 + e * u1 + e * e * u2
    v = math.sqrt(0.5 * (gamma + 1.0) * e) * (e * v1 + e * e * v2)
    return u, v


def hall_line(
    gamma: float, Rt: float, Rc: float, *, axis_mach: float, n: int, delta: float = 1.0
) -> tuple[float, float, list[tuple[float, float, float]]]:
    """Hall 2nd-order initial-data line; same contract as :func:`sauer_line`.

    Axisymmetric only (``ValueError`` for planar, so callers fall back to Sauer).
    """
    if abs(delta) < 0.5:
        raise ValueError("Hall series implemented for axisymmetric throats only")
    if Rc / Rt < HALL_MIN_RC_RT:
        raise ValueError(f"Hall series diverges for Rc/Rt < {HALL_MIN_RC_RT}; use Kliegel–Levine or Sauer")
    return _series_line(hall_uv, gamma, Rt, Rc, axis_mach=axis_mach, n=n)


def kl_line(
    gamma: float, Rt: float, Rc: float, *, axis_mach: float, n: int, delta: float = 1.0
) -> tuple[float, float, list[tuple[float, float, float]]]:
    """Kliegel–Levine initial-data line; same contract as :func:`sauer_line`."""
    if abs(delta) < 0.5:
        raise ValueError("Kliegel–Levine series implemented for axisymmetric throats only")
    if Rc / Rt < KL_MIN_RC_RT:
        raise ValueError(f"Kliegel–Levine series unvalidated for Rc/Rt < {KL_MIN_RC_RT}; use Sauer")
    return _series_line(kl_uv, gamma, Rt, Rc, axis_mach=axis_mach, n=n)


def _series_line(uv, gamma: float, Rt: float, Rc: float, *, axis_mach: float, n: int):
    """Vertical start line where the axis reaches ``axis_mach`` for a (u, v) series."""
    target = mstar_from_mach(axis_mach, gamma)
    try:
        x = sauer_line(gamma, Rt, Rc, axis_mach=axis_mach, n=1)[0]
    except ValueError:  # Sauer's line can lie past the arc for tight throats
        x = 0.1 * Rt
    for _ in range(50):
        f = uv(x, 0.0, gamma, Rt, Rc)[0] - target
        h = 1e-7 * Rt
        df = (uv(x + h, 0.0, gamma, Rt, Rc)[0] - f - target) / h
        if df <= 0.0:
            raise ValueError("series axis velocity not increasing (Rc/Rt too small)")
        x -= f / df
        if abs(f) < 1e-13:
            break
    else:
        raise ValueError("transonic initial line did not converge")
    if not 0.0 < x < Rc:
        raise ValueError("axis_mach too high for this throat radius of curvature")
    phi = math.asin(x / Rc)
    r_w = Rt + Rc * (1.0 - math.cos(phi))
    pts = []
    for j in range(n + 1):
        r = r_w * j / n
        u, v = uv(x, r, gamma, Rt, Rc)
        th = phi if j == n else math.atan2(v, u)
        pts.append((r, mach_from_mstar(math.hypot(u, v), gamma), th))
    return x, phi, pts


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
