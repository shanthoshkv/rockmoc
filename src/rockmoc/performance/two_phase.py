"""Gas–particle (two-phase) nozzle loss: quasi-1-D velocity and thermal lag.

Condensed particles (e.g. Al2O3, mass ratio κ = ṁ_p/ṁ_g, diameter d) exchange
momentum and heat with the gas along the prescribed area A(x)::

    du_p/dx = f_D (u − u_p) / (τ_v u_p),   τ_v = ρ_p d² / (18 μ),   f_D = 1 + 0.15 Re_p^0.687
    dT_p/dx = Nu/2 · (T − T_p) / (τ_T u_p), τ_T = ρ_p c_s d² / (12 k), Nu = 2 + 0.6 Re_p^½ Pr^⅓

and the gas obeys continuity ρuA = ṁ_g, mixture momentum and energy, giving

    du/dx = γ / (1 − 1/M²) · [ (RT/u) A'/A + (R/(c_p u)) E_p − D_p ]
    D_p = κ du_p/dx,  E_p = κ (c_s dT_p/dx + u_p du_p/dx),  dT/dx = −(u du/dx + E_p)/c_p

(derived from conservation; Hoglund 1962 framework). The no-lag reference is
the equilibrium mixture isentrope (γ_m, R_m); η_2φ = Isp_lag / Isp_no-lag.
Integration starts in equilibrium at a supersonic area ratio (default 1.05)
— ponytail: real particles already lag in the convergent section, so η_2φ is
slightly optimistic; a full subsonic–supersonic solve would remove that.
Particle size is an input: the Hermsen D43 correlation was not verified.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from rockmoc.thermo.perfect_gas import PerfectGas

G0 = 9.80665


@dataclass
class TwoPhaseResult:
    Isp_lag: float
    Isp_no_lag: float
    eta_two_phase: float  # pass as eta_two_phase
    x: np.ndarray
    u: np.ndarray
    u_p: np.ndarray
    T: np.ndarray
    T_p: np.ndarray


def two_phase_nozzle(
    *,
    x: np.ndarray,
    r: np.ndarray,
    gamma_gas: float,
    R_gas: float,
    pc: float,
    Tc: float,
    particle_mass_fraction: float,
    diameter: float,
    rho_p: float = 3970.0,
    c_s: float = 1380.0,
    mu_gas: float = 9.0e-5,
    prandtl: float = 0.8,
    start_area_ratio: float = 1.05,
) -> TwoPhaseResult:
    """Two-phase expansion along r(x); ``particle_mass_fraction`` = ṁ_p/(ṁ_g + ṁ_p)."""
    if not 0.0 <= particle_mass_fraction < 1.0:
        raise ValueError("particle_mass_fraction must be in [0, 1)")
    kap = particle_mass_fraction / (1.0 - particle_mass_fraction)
    g, R = gamma_gas, R_gas
    cp = g * R / (g - 1.0)
    k_gas = mu_gas * cp / prandtl
    # equilibrium (no-lag) mixture
    cpm = (cp + kap * c_s) / (1.0 + kap)
    Rm = R / (1.0 + kap)
    mix = PerfectGas(cpm / (cpm - Rm), Rm)

    x = np.asarray(x, dtype=np.float64)
    r = np.asarray(r, dtype=np.float64)
    Rt = float(np.min(r))
    eps_x = (r / Rt) ** 2
    i_t = int(np.argmin(r))
    i0 = int(np.argmax((np.arange(len(x)) >= i_t) & (eps_x >= start_area_ratio)))
    xs, A = x[i0:], math.pi * r[i0:] ** 2
    dAdx = np.gradient(A, xs)

    def mix_state(eps: float) -> tuple[float, float, float]:
        M = float(mix.mach_from_area_ratio(eps, branch="supersonic"))
        T = Tc * float(1.0 / mix.temperature_ratio(M))
        p = pc * float(mix.pe_pc_from_mach(M))
        return T, p, M * math.sqrt(mix.gamma * Rm * T)

    T0_, p0, u0 = mix_state(float(eps_x[i0]))
    mdot_g = p0 / (R * T0_) * u0 * A[0]
    mdot = (1.0 + kap) * mdot_g

    def rhs(xi: float, y: np.ndarray) -> np.ndarray:
        u, T, up, Tp = y
        Ai = float(np.interp(xi, xs, A))
        rho = mdot_g / (u * Ai)
        Re = max(rho * abs(u - up) * diameter / mu_gas, 0.0)
        tv = rho_p * diameter**2 / (18.0 * mu_gas)
        tT = rho_p * c_s * diameter**2 / (12.0 * k_gas)
        dup = (1.0 + 0.15 * Re**0.687) * (u - up) / (tv * up) if kap > 0 else 0.0
        dTp = (2.0 + 0.6 * math.sqrt(Re) * prandtl ** (1 / 3)) / 2.0 * (T - Tp) / (tT * up) if kap > 0 else 0.0
        Dp = kap * dup
        Ep = kap * (c_s * dTp + up * dup)
        M2 = u * u / (g * R * T)
        du = g / (1.0 - 1.0 / M2) * (R * T / u * float(np.interp(xi, xs, dAdx)) / Ai + R / (cp * u) * Ep - Dp)
        dT = -(u * du + Ep) / cp
        return np.array([du, dT, dup, dTp])

    sol = solve_ivp(rhs, (float(xs[0]), float(xs[-1])), [u0, T0_, u0, T0_], method="LSODA",
                    t_eval=xs, rtol=1e-8, atol=1e-8)
    if not sol.success:
        raise RuntimeError(f"two-phase integration failed: {sol.message}")
    u, T, up, Tp = sol.y
    pe = mdot_g / (u[-1] * A[-1]) * R * T[-1]
    isp_lag = ((mdot_g * u[-1] + kap * mdot_g * up[-1]) + pe * A[-1]) / (mdot * G0)
    Te, pe0, ue = mix_state(float(eps_x[-1]))
    isp_eq = (ue + pe0 * A[-1] / mdot) / G0
    return TwoPhaseResult(isp_lag, isp_eq, isp_lag / isp_eq, xs, u, up, T, Tp)
