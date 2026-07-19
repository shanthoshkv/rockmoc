"""Characteristic slopes and axisymmetric compatibility — first principles.

Theory sources (local copies under ``refs/``)
--------------------------------------------
* Zucrow & Hoffman, *Gas Dynamics* — axisymmetric MOC unit processes
* Shapiro, *Dynamics and Thermodynamics of Compressible Fluid Flow*
* NASA TM-X-1502 — planar sharp-throat MLN (averaged slopes)
* Young (2012) / classical nozzle MOC notes — axisym source terms

Irrotational isentropic perfect-gas flow
---------------------------------------
Prandtl–Meyer differential (along an isentrope)::

    dν = √(M² − 1)  dV / V

Characteristic directions (meridional plane)::

    (dr/dx)_± = tan(θ ± μ),   μ = arcsin(1/M)

**Planar (δ = 0)** Riemann invariants::

    along C+:  θ + ν = constant
    along C−:  θ − ν = constant

**Axisymmetric (δ = 1)** geometric source (Zucrow/Hoffman / standard form)::

    along C+ [slope tan(θ+μ)]:
        d(θ + ν)/dx = − δ · sinθ · sinμ / (r · cos(θ+μ))

    along C− [slope tan(θ−μ)]:
        d(θ − ν)/dx = − δ · sinθ · sinμ / (r · cos(θ−μ))

Finite-difference unit processes use **averaged** angles/radii between the
known point and the unknown point (predictor–corrector), as in TM-X-1502
and Hoffman-style codes.

Axis limit
----------
As r → 0, θ → 0 the geometric source is singular in appearance but the
product sinθ/r remains finite for smooth solutions. Implementations use a
small-r cutoff and special axis unit processes (θ_P = 0).
"""

from __future__ import annotations

import math

from rockmoc.moc.types import FlowPoint

_R_EPS = 1e-14
_ANGLE_EPS = 1e-14


def char_slope_plus(theta: float, mu: float) -> float:
    """dr/dx along C+ characteristic."""
    return math.tan(theta + mu)


def char_slope_minus(theta: float, mu: float) -> float:
    """dr/dx along C− characteristic."""
    return math.tan(theta - mu)


def _source_core(theta: float, mu: float, r: float, char_angle: float, delta: float) -> float:
    """Shared axisymmetric source −δ sinθ sinμ / (r cos φ), φ = θ±μ.

    Near the axis, θ ∼ O(r) for smooth solutions so sinθ/r stays finite.
    We regularize with r_eff = max(r, r_floor) and sinθ → θ when |θ| is small.
    """
    if abs(delta) < _ANGLE_EPS:
        return 0.0
    r_floor = 1e-10
    r_eff = r if abs(r) > r_floor else r_floor
    # sinθ / r : use θ/r when |θ| small (consistent with axis θ→0)
    if abs(theta) < 1e-8:
        sin_theta_over_r = theta / r_eff  # ≈ 0 on axis if θ=0
    else:
        sin_theta_over_r = math.sin(theta) / r_eff
    cos_phi = math.cos(char_angle)
    if abs(cos_phi) < 1e-12:
        return 0.0
    return -float(delta) * sin_theta_over_r * math.sin(mu) / cos_phi


def source_plus(theta: float, mu: float, r: float, *, delta: float = 1.0) -> float:
    """Geometric source S+ so that (θ+ν)_P ≈ (θ+ν)_B + S+ · (x_P − x_B).

    First-principles axisymmetric form (δ = 1)::

        S+ = − sinθ · sinμ / (r · cos(θ+μ))
    """
    return _source_core(theta, mu, r, theta + mu, delta)


def source_minus(theta: float, mu: float, r: float, *, delta: float = 1.0) -> float:
    """Geometric source S− so that (θ−ν)_P ≈ (θ−ν)_A + S− · (x_P − x_A).

    First-principles::

        S− = − sinθ · sinμ / (r · cos(θ−μ))
    """
    return _source_core(theta, mu, r, theta - mu, delta)

def source_along_ds(
    theta: float,
    mu: float,
    r: float,
    *,
    family: str,
    delta: float = 1.0,
) -> float:
    """Source for d(θ±ν)/ds where s is arc length along the characteristic.

    ds² = dx² + dr²,  dr/dx = tan(θ±μ) ⇒ ds/dx = |sec(θ±μ)|.
    So d(·)/dx = d(·)/ds · ds/dx.

    Returns d(θ±ν)/ds coefficient for use when marching in arc length.
    """
    if abs(delta) < _ANGLE_EPS or abs(r) < _R_EPS:
        return 0.0
    if family == "+":
        # d(θ+ν)/dx = S+ ⇒ d(θ+ν)/ds = S+ / (ds/dx)
        ang = theta + mu
        Sx = source_plus(theta, mu, r, delta=delta)
    else:
        ang = theta - mu
        Sx = source_minus(theta, mu, r, delta=delta)
    ds_dx = abs(1.0 / max(abs(math.cos(ang)), 1e-12))
    return Sx / ds_dx


def average_state(a: FlowPoint, b: FlowPoint) -> tuple[float, float, float, float]:
    """Average θ, μ, r, ν for corrector steps (Hoffman / TM-X-1502 style)."""
    return (
        0.5 * (a.theta + b.theta),
        0.5 * (a.mu + b.mu),
        0.5 * (a.r + b.r),
        0.5 * (a.nu + b.nu),
    )


def riemann_invariants(theta: float, nu: float) -> tuple[float, float]:
    """Return (J+, J−) = (θ+ν, θ−ν). Planar invariants when δ=0."""
    return theta + nu, theta - nu


def from_riemann(jp: float, jm: float) -> tuple[float, float]:
    """θ, ν from Riemann invariants J±."""
    theta = 0.5 * (jp + jm)
    nu = 0.5 * (jp - jm)
    return theta, nu
