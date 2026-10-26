"""Quasi-1-D finite-rate (kinetic) nozzle expansion with Cantera — ODK-style.

Chamber: HP equilibrium at p_c. Shifting equilibrium through the throat up to a
supersonic start area ratio (default 1.05, avoids the frozen-sonic singularity),
then finite-rate chemistry along the wall area A(x)::

    dY_k/dx = ω̇_k W_k / (ρ u)
    du/dx   = u / (1 − M_f²) · [ −A'/A + Σ_k (W̄/W_k − h_k/(c_p T)) dY_k/dx ]
    dT/dx   = −( u du/dx + Σ_k h_k dY_k/dx ) / c_p          (h_k per unit mass)

with ρ = ṁ/(u A) and M_f the frozen Mach number. Derived from conservation
(the ODK/TDK reaction sets themselves were not available); checked in the
tests by its limits: zero rates reproduce the frozen expansion and fast rates
the shifting-equilibrium one. Kinetic loss matters for small, low-p_c, high-ε
and hydrogen engines (RocketIsp: 0.1–10 % of Isp; throat size dominates).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

G0 = 9.80665


def _require_cantera():
    try:
        import cantera as ct  # type: ignore
    except ImportError as exc:
        raise ImportError("Cantera is required: pip install rockmoc[cantera]") from exc
    return ct


@dataclass
class KineticResult:
    Isp_kinetic: float  # vacuum, s
    Isp_shifting: float
    Isp_frozen: float  # composition frozen at the kinetic start station
    eta_kinetic: float  # Isp_kinetic / Isp_shifting → pass as eta_chemical
    frac_kinetic: float  # (kin − frozen)/(shifting − frozen), RocketIsp definition
    x: np.ndarray
    T: np.ndarray
    p: np.ndarray
    u: np.ndarray
    X_exit: dict


def kinetic_nozzle(
    *,
    mechanism: str = "h2o2.yaml",
    reactants: str = "H2:2.0, O2:1.0",
    T_reactants: float = 300.0,
    pc: float = 7e6,
    x: np.ndarray,
    r: np.ndarray,
    start_area_ratio: float = 1.05,
    rate_multiplier: float = 1.0,
) -> KineticResult:
    """Finite-rate expansion along a wall r(x) (x from the throat, r[0] at or near the throat)."""
    ct = _require_cantera()
    gas = ct.Solution(mechanism)
    gas.TPX = T_reactants, pc, reactants
    gas.equilibrate("HP")
    h0, s0, Yc = gas.enthalpy_mass, gas.entropy_mass, gas.Y.copy()

    def eq_state(p: float, frozen_Y: np.ndarray | None = None) -> tuple[float, float, float]:
        if frozen_Y is None:
            gas.SPY = s0, p, Yc
            gas.equilibrate("SP")
        else:
            gas.SPY = s_f, p, frozen_Y
        u = math.sqrt(max(2.0 * (h_ref - gas.enthalpy_mass), 0.0))
        return gas.density, u, gas.enthalpy_mass

    h_ref = h0
    s_f = s0
    best = minimize_scalar(lambda p: -eq_state(p)[0] * eq_state(p)[1], bounds=(0.3 * pc, 0.9 * pc),
                           method="bounded", options={"xatol": 1e-6 * pc})
    p_star = float(best.x)
    rho_s, u_s, _ = eq_state(p_star)
    G_star = rho_s * u_s

    def p_at(eps: float, frozen_Y: np.ndarray | None = None) -> float:
        def f(p: float) -> float:
            rho, u, _ = eq_state(p, frozen_Y)
            return G_star / max(rho * u, 1e-300) - eps
        return float(brentq(f, 1e-7 * pc, p_star * 0.999, xtol=1e-10 * pc))

    def isp(rho: float, u: float, p: float) -> float:
        return (u + p / (rho * u)) / G0

    x = np.asarray(x, dtype=np.float64)
    r = np.asarray(r, dtype=np.float64)
    Rt = float(np.min(r))
    eps_x = (r / Rt) ** 2
    eps_e = float(eps_x[-1])
    i_t = int(np.argmin(r))
    # start station on the supersonic side
    sup = np.arange(len(x)) >= i_t
    i0 = int(np.argmax(sup & (eps_x >= start_area_ratio)))
    x0 = float(x[i0])

    p0 = p_at(float(eps_x[i0]))
    rho0, u0, _ = eq_state(p0)
    Y0, T0_ = gas.Y.copy(), gas.T
    # shifting equilibrium and frozen-from-start references at the exit
    pe = p_at(eps_e)
    rho_e, u_e, _ = eq_state(pe)
    isp_eq = isp(rho_e, u_e, pe)
    gas.TPY = T0_, p0, Y0
    s_f = gas.entropy_mass
    pf = p_at(eps_e, frozen_Y=Y0)
    rho_f, u_f, _ = eq_state(pf, frozen_Y=Y0)
    isp_fr = isp(rho_f, u_f, pf)

    mdot = G_star * math.pi * Rt * Rt
    xs, rs = x[i0:], r[i0:]
    A = math.pi * rs**2
    dAdx = np.gradient(A, xs)
    W = gas.molecular_weights
    gas.set_multiplier(rate_multiplier)

    def rhs(xi: float, y: np.ndarray) -> np.ndarray:
        u, T = y[0], y[1]
        Y = np.clip(y[2:], 0.0, None)
        Ai = float(np.interp(xi, xs, A))
        gas.TDY = T, mdot / (u * Ai), Y
        rho, cp = gas.density, gas.cp_mass
        dY = gas.net_production_rates * W / (rho * u)
        hk = gas.partial_molar_enthalpies / W
        a2 = gas.cp_mass / gas.cv_mass * ct.gas_constant / gas.mean_molecular_weight * T
        Mf2 = u * u / a2
        S = float(np.sum((gas.mean_molecular_weight / W - hk / (cp * T)) * dY))
        du = u / (1.0 - Mf2) * (-float(np.interp(xi, xs, dAdx)) / Ai + S)
        dT = -(u * du + float(np.sum(hk * dY))) / cp
        return np.concatenate([[du, dT], dY])

    sol = solve_ivp(rhs, (x0, float(xs[-1])), np.concatenate([[u0, T0_], Y0]), method="BDF",
                    t_eval=xs, rtol=1e-7, atol=1e-10)
    gas.set_multiplier(1.0)
    if not sol.success:
        raise RuntimeError(f"kinetic integration failed: {sol.message}")
    u, T = sol.y[0], sol.y[1]
    Ye = np.clip(sol.y[2:, -1], 0.0, None)
    rho_k = mdot / (u[-1] * A[-1])
    gas.TDY = T[-1], rho_k, Ye
    isp_k = isp(rho_k, u[-1], gas.P)
    p = np.empty(len(u))
    for i, (Ti, ui, Ai, Yi) in enumerate(zip(T, u, A, sol.y[2:].T)):
        gas.TDY = Ti, mdot / (ui * Ai), np.clip(Yi, 0.0, None)
        p[i] = gas.P
    gas.TDY = T[-1], rho_k, Ye
    frac = (isp_k - isp_fr) / (isp_eq - isp_fr) if isp_eq > isp_fr else 1.0
    return KineticResult(isp_k, isp_eq, isp_fr, isp_k / isp_eq, frac, xs, T, p, u,
                         {n: float(v) for n, v in zip(gas.species_names, gas.X) if v > 1e-4})
