"""Isentropic expansion property tables (variable / multi-γ support).

Builds a table along s = const from chamber stagnation to high expansion
using either perfect-gas closed form or optional CEA equilibrium sampling.

MOC can query local a, ρ, γ_eff, T, p vs Mach or pressure ratio.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from rockmoc.exceptions import ThermoError
from rockmoc.thermo.perfect_gas import PerfectGas


@dataclass
class ExpansionTable:
    """Tabulated isentropic expansion states."""

    gamma_ref: float
    R: float
    pc: float
    Tc: float
    M: NDArray[np.floating]
    p_pc: NDArray[np.floating]
    T_Tc: NDArray[np.floating]
    rho_rhoc: NDArray[np.floating]
    a: NDArray[np.floating]  # m/s
    gamma_local: NDArray[np.floating]
    nu: NDArray[np.floating]  # rad, integrated generalized PM
    area_ratio: NDArray[np.floating]
    mode: str
    meta: dict

    def _interp(self, name: str, M: float) -> float:
        y = getattr(self, name)
        return float(np.interp(M, self.M, y))

    def state_at_mach(self, M: float) -> dict[str, float]:
        M = float(np.clip(M, self.M[0], self.M[-1]))
        return {
            "M": M,
            "p_pc": self._interp("p_pc", M),
            "T": self.Tc * self._interp("T_Tc", M),
            "rho": (self.pc / (self.R * self.Tc)) * self._interp("rho_rhoc", M),
            "a": self._interp("a", M),
            "gamma": self._interp("gamma_local", M),
            "nu": self._interp("nu", M),
            "area_ratio": self._interp("area_ratio", M),
        }

    def mach_from_nu(self, nu: float) -> float:
        nu = float(np.clip(nu, self.nu[0], self.nu[-1]))
        return float(np.interp(nu, self.nu, self.M))

    def mach_from_area_ratio(self, eps: float) -> float:
        """Supersonic branch inverse A/A*."""
        eps = max(float(eps), 1.0)
        # find on supersonic part M>=1
        return float(np.interp(eps, self.area_ratio, self.M))

    def as_perfect_gas(self) -> PerfectGas:
        """Frozen reference gas using table γ_ref and R."""
        return PerfectGas(gamma=self.gamma_ref, R=self.R)


def build_perfect_gas_table(
    gas: PerfectGas,
    *,
    pc: float = 5e6,
    Tc: float = 3000.0,
    M_max: float = 8.0,
    n: int = 400,
) -> ExpansionTable:
    """Isentropic table for constant-γ perfect gas with integrated ν(M)."""
    if M_max <= 1.0 or n < 20:
        raise ThermoError("M_max > 1 and n ≥ 20 required")
    M = np.linspace(1.0, M_max, n)
    g = gas.gamma
    T_ratio = 1.0 / (1.0 + 0.5 * (g - 1.0) * M**2)
    p_ratio = T_ratio ** (g / (g - 1.0))
    rho_ratio = T_ratio ** (1.0 / (g - 1.0))
    T = Tc * T_ratio
    a = np.sqrt(g * gas.R * T)
    nu = np.array([float(gas.prandtl_meyer(m)) for m in M])
    ar = np.array([float(gas.area_ratio(m)) for m in M])
    gamma_loc = np.full_like(M, g)
    return ExpansionTable(
        gamma_ref=g,
        R=gas.R,
        pc=pc,
        Tc=Tc,
        M=M,
        p_pc=p_ratio,
        T_Tc=T_ratio,
        rho_rhoc=rho_ratio,
        a=a,
        gamma_local=gamma_loc,
        nu=nu,
        area_ratio=ar,
        mode="perfect_gas",
        meta={"n": n, "M_max": M_max},
    )


def build_multi_gamma_table(
    *,
    gamma_c: float,
    gamma_t: float,
    gamma_e: float,
    R: float,
    pc: float = 5e6,
    Tc: float = 3000.0,
    M_max: float = 8.0,
    n: int = 400,
) -> ExpansionTable:
    """Thermally perfect isentrope with γ(M) blended chamber → throat → exit.

    γ is linear in M from γ_c (M = 0) to γ_t (M = 1) and from γ_t to γ_e
    (M = M_max). The expansion is integrated in T from stagnation with
    c_p = γR/(γ−1) (frozen composition, calorically imperfect)::

        dh = c_p dT,   d ln p = c_p dT / (R T),   V² = 2 (h_c − h)
        a² = γ R T,    A/A* = (ρV)* / (ρV),       dν = √(M² − 1) dV / V

    (ρV peaks exactly at M = 1 for any such gas). This is consistent where the
    old local-γ closed forms were not; for γ_c = γ_t = γ_e it reproduces the
    perfect-gas relations.
    """
    if min(gamma_c, gamma_t, gamma_e) <= 1.0:
        raise ThermoError("all gamma values must be > 1")
    if M_max <= 1.0 or n < 20:
        raise ThermoError("M_max > 1 and n ≥ 20 required")

    def g_of(m: float) -> float:
        if m <= 1.0:
            return gamma_c + (gamma_t - gamma_c) * m
        return gamma_t + (gamma_e - gamma_t) * min((m - 1.0) / (M_max - 1.0), 1.0)

    # Integrate in t = T/Tc (R = 1 per unit Tc; ratios do not depend on R, Tc)
    g_min = min(gamma_c, gamma_t, gamma_e)
    t_end = 0.5 / (1.0 + 0.5 * (g_min - 1.0) * M_max**2)
    n_int = 40 * n
    ts = np.geomspace(1.0, t_end, n_int)
    h = np.zeros(n_int)  # (h_c − h) / (R Tc)
    lnp = np.zeros(n_int)
    Ms = np.zeros(n_int)
    gs = np.full(n_int, gamma_c)
    for i in range(1, n_int):
        dt = ts[i] - ts[i - 1]
        cp0 = gs[i - 1] / (gs[i - 1] - 1.0)
        g1 = gs[i - 1]
        for _ in range(2):  # predictor–corrector on γ at the new station
            cp = 0.5 * (cp0 + g1 / (g1 - 1.0))
            h_i = h[i - 1] - cp * dt
            m_i = math.sqrt(2.0 * h_i / (g1 * ts[i]))
            g1 = g_of(m_i)
        h[i], Ms[i], gs[i] = h_i, m_i, g1
        lnp[i] = lnp[i - 1] + cp * math.log(ts[i] / ts[i - 1])
        if m_i > M_max:
            break
    k = i + 1
    ts, h, lnp, Ms, gs = ts[:k], h[:k], lnp[:k], Ms[:k], gs[:k]
    if Ms[-1] < M_max:
        raise ThermoError("isentrope integration did not reach M_max")
    p_r = np.exp(lnp)
    V = np.sqrt(2.0 * h)
    G = p_r / ts * V  # ρV / (ρ_c √(R Tc))

    M = np.linspace(1.0, M_max, n)
    T_ratio = np.interp(M, Ms, ts)
    p_ratio = np.exp(np.interp(M, Ms, lnp))
    rho_ratio = p_ratio / T_ratio
    gamma = np.array([g_of(float(m)) for m in M])
    a = np.sqrt(gamma * R * Tc * T_ratio)
    G_star = float(np.interp(1.0, Ms, G))
    ar = G_star / np.interp(M, Ms, G)
    ar[0] = 1.0
    sup = Ms >= 1.0
    Vs, Msup = V[sup], Ms[sup]
    dnu = np.sqrt(np.maximum(Msup**2 - 1.0, 0.0)) / Vs
    # first supersonic sample sits within ~1e-4 of M = 1, so ν(it) ≈ 0
    nu_s = np.concatenate([[0.0], np.cumsum(0.5 * (dnu[1:] + dnu[:-1]) * np.diff(Vs))])
    nu = np.maximum(np.interp(M, Msup, nu_s), 0.0)
    nu[0] = 0.0
    # c* = p_c A* / ṁ = p_c / (ρV)* in absolute units
    cstar = float(math.sqrt(R * Tc) / G_star)

    return ExpansionTable(
        gamma_ref=gamma_t,
        R=R,
        pc=pc,
        Tc=Tc,
        M=M,
        p_pc=p_ratio,
        T_Tc=T_ratio,
        rho_rhoc=rho_ratio,
        a=a,
        gamma_local=gamma,
        nu=nu,
        area_ratio=ar,
        mode="multi_gamma",
        meta={
            "cstar_m_s": cstar,
            "gamma_c": gamma_c,
            "gamma_t": gamma_t,
            "gamma_e": gamma_e,
            "n": n,
            "M_max": M_max,
        },
    )


def build_cea_frozen_table(
    ox: str,
    fuel: str,
    of_ratio: float,
    pc_psia: float,
    *,
    eps_max: float = 100.0,
    n: int = 50,
) -> ExpansionTable:
    """Build expansion table from RocketCEA frozen performance samples.

    Samples several expansion ratios with frozen composition and interpolates.
    Requires ``pip install rockmoc[cea]``.
    """
    try:
        from rocketcea.cea_obj import CEA_Obj  # type: ignore
    except ImportError as exc:
        raise ImportError("RocketCEA required: pip install rockmoc[cea]") from exc

    cea = CEA_Obj(oxName=ox, fuelName=fuel)
    Pc = float(pc_psia)
    ofr = float(of_ratio)
    eps_list = np.geomspace(1.01, eps_max, n)

    M_list, p_list, T_list, g_list, ar_list = [], [], [], [], []
    Tc = float(cea.get_Tcomb(Pc, ofr))
    mw, gamma_c = cea.get_Chamber_MolWt_gamma(Pc, ofr, float(eps_max))
    try:
        mw_t, gamma_t = cea.get_Throat_MolWt_gamma(Pc, ofr, float(eps_max))
    except Exception:
        mw_t, gamma_t = mw, gamma_c
    R = 8314.462618 / float(mw_t)

    for eps in eps_list:
        try:
            # frozen Isp / Mach estimates
            Me = float(cea.get_MachNumber(Pc, ofr, float(eps), frozen=1, frozenAtThroat=1))
            # temperatures
            # RocketCEA: get_Temperatures returns chamber, throat, exit
            temps = cea.get_Temperatures(Pc, ofr, float(eps), frozen=1, frozenAtThroat=1)
            Te = float(temps[-1]) if temps else Tc / (1 + 0.2 * Me**2)
            # pressures via isentropic with local gamma
            g_use = float(gamma_t)
            pe_pc = (1.0 + 0.5 * (g_use - 1.0) * Me**2) ** (-g_use / (g_use - 1.0))
        except Exception:
            gas = PerfectGas(gamma=float(gamma_t), R=R)
            Me = float(gas.mach_from_area_ratio(float(eps), branch="supersonic"))
            pe_pc = float(gas.pe_pc_from_mach(Me))
            Te = Tc / float(gas.temperature_ratio(Me))
            g_use = float(gamma_t)
        M_list.append(max(Me, 1.0))
        p_list.append(pe_pc)
        T_list.append(Te / Tc)
        g_list.append(g_use)
        ar_list.append(float(eps))

    # Sort by M
    order = np.argsort(M_list)
    M = np.array(M_list)[order]
    p_pc = np.array(p_list)[order]
    T_Tc = np.array(T_list)[order]
    gamma = np.array(g_list)[order]
    ar = np.array(ar_list)[order]
    # densify unique M
    _, uniq = np.unique(np.round(M, 6), return_index=True)
    M, p_pc, T_Tc, gamma, ar = M[uniq], p_pc[uniq], T_Tc[uniq], gamma[uniq], ar[uniq]
    T = Tc * T_Tc
    a = np.sqrt(gamma * R * T)
    rho_rhoc = p_pc / T_Tc  # p/ρ ~ T for ideal
    gas = PerfectGas(gamma=float(gamma_t), R=R)
    nu = np.array([float(gas.prandtl_meyer(m)) if m >= 1 else 0.0 for m in M])

    return ExpansionTable(
        gamma_ref=float(gamma_t),
        R=R,
        pc=pc_psia * 6894.757,  # store Pa
        Tc=Tc,
        M=M,
        p_pc=p_pc,
        T_Tc=T_Tc,
        rho_rhoc=rho_rhoc,
        a=a,
        gamma_local=gamma,
        nu=nu,
        area_ratio=ar,
        mode="cea_frozen",
        meta={
            "ox": ox,
            "fuel": fuel,
            "OF": ofr,
            "Pc_psia": Pc,
            "gamma_c": float(gamma_c),
            "gamma_t": float(gamma_t),
            "mw": float(mw_t),
        },
    )
