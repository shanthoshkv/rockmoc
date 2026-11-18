"""Real-gas isentropic expansion table from RocketCEA (equilibrium or frozen-at-throat).

Unlike a constant-γ gas, every state along the nozzle comes from the thermochemistry
solver: for each area ratio ε CEA returns p, T, ρ, a, h and M. From those::

    γ_s = ρ a² / p                       local isentropic exponent (exact for ρV² = γ_s p M²)
    V   = √(2 (h_c − h))                 energy equation (checked against M·a, see meta)
    ν   = ∫ √(M² − 1) d ln V             generalized Prandtl–Meyer angle (any gas)
    A/A* = ε                             by definition of CEA's throat
    c*  = p_c / (ρ a)*                   at M = 1

so the MOC unit processes, the Prandtl–Meyer start, pe/pc ↔ M ↔ A/A* and the Cf/Isp
formulas of :class:`~rockmoc.thermo.variable_gamma.VariableGammaGas` all see the same
variable-γ, variable-molecular-weight gas. ``mode``:

* ``"equilibrium"``  shifting equilibrium (upper bound on Isp; ideal fast chemistry)
* ``"frozen"``       composition frozen at the throat (lower bound; no recombination)

Real engines lie between; sweep ``gas_model=frozen,equilibrium`` to bracket them.
Requires ``pip install rocketcea``.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.interpolate import PchipInterpolator

from rockmoc.exceptions import ThermoError
from rockmoc.thermo.expansion_table import ExpansionTable

PSI = 6894.757293168  # Pa
_FT = 0.3048  # ft → m
_LBFT3 = 16.01846337  # lbm/ft³ → kg/m³
_BTULB = 2326.0  # BTU/lbm → J/kg
MODES = ("equilibrium", "frozen")


def build_cea_table(
    ox: str,
    fuel: str,
    of: float,
    pc: float,
    *,
    mode: str = "equilibrium",
    eps_max: float = 400.0,
    n: int = 70,
) -> ExpansionTable:
    """Expansion table for propellants ``ox``/``fuel`` at O/F ``of`` and chamber ``pc`` [Pa].

    ``eps_max`` bounds the table (area ratio); designs beyond it are rejected, not clipped.
    Results are cached per argument set.
    """
    if mode not in MODES:
        raise ThermoError(f"mode must be one of {MODES}; got {mode!r}")
    if eps_max <= 2.0 or n < 20:
        raise ThermoError("eps_max > 2 and n ≥ 20 required")
    return _build(ox, fuel, float(of), float(pc), mode, float(eps_max), int(n))


@lru_cache(maxsize=64)
def _build(ox: str, fuel: str, of: float, pc: float, mode: str, eps_max: float, n: int):
    try:
        from rocketcea.cea_obj import CEA_Obj
    except ImportError as exc:
        raise ImportError("RocketCEA required: pip install rocketcea") from exc

    cea = CEA_Obj(oxName=ox, fuelName=fuel)
    Pc = pc / PSI  # psia
    fz = 1 if mode == "frozen" else 0
    kw = {"frozen": fz, "frozenAtThroat": fz}

    eps = np.geomspace(1.003, eps_max, n)
    M, pr, T, a, rho, h = ([] for _ in range(6))
    chamber = throat = None
    for e in eps:
        e = float(e)
        Ms = float(cea.get_MachNumber(Pc, of, e, **kw))
        pcpe = float(cea.get_PcOvPe(Pc, of, e, **kw))
        Ts = cea.get_Temperatures(Pc, of, e, **kw)
        as_ = cea.get_SonicVelocities(Pc, of, e, **kw)
        rs = cea.get_Densities(Pc, of, e, **kw)
        hs = cea.get_Enthalpies(Pc, of, e, **kw)
        if chamber is None:
            chamber = (Ts[0] / 1.8, as_[0] * _FT, rs[0] * _LBFT3, hs[0] * _BTULB)
            throat = (Ts[1] / 1.8, as_[1] * _FT, rs[1] * _LBFT3, hs[1] * _BTULB)
        M.append(Ms)
        pr.append(1.0 / pcpe)
        T.append(Ts[2] / 1.8)
        a.append(as_[2] * _FT)
        rho.append(rs[2] * _LBFT3)
        h.append(hs[2] * _BTULB)
    assert chamber is not None and throat is not None
    Tc, ac, rhoc, hc = chamber
    Tt, at, rhot, ht = throat
    pt_pc = 1.0 / float(cea.get_Throat_PcOvPe(Pc, of))  # same throat in both modes

    # prepend the throat (M = 1, ε = 1)
    M_a = np.concatenate([[1.0], M])
    p_pc = np.concatenate([[pt_pc], pr])
    T_a = np.concatenate([[Tt], T])
    a_a = np.concatenate([[at], a])
    rho_a = np.concatenate([[rhot], rho])
    h_a = np.concatenate([[ht], h])
    eps_a = np.concatenate([[1.0], eps])
    keep = np.concatenate([[True], np.diff(M_a) > 1e-9])  # strictly increasing M
    M_a, p_pc, T_a, a_a, rho_a, h_a, eps_a = (
        x[keep] for x in (M_a, p_pc, T_a, a_a, rho_a, h_a, eps_a)
    )

    p_abs = p_pc * pc
    gamma = rho_a * a_a**2 / p_abs
    V = np.sqrt(np.maximum(2.0 * (hc - h_a), 0.0))
    V[0] = at
    resid = float(np.max(np.abs(V[1:] / (M_a[1:] * a_a[1:]) - 1.0)))

    # ν(M) = ∫ √(M²−1) d ln V — smooth PCHIP of ln V(M), fine quadrature
    lnV = PchipInterpolator(M_a, np.log(V))
    dlnV = lnV.derivative()
    Md = np.linspace(1.0, M_a[-1], 20000)
    f = np.sqrt(np.maximum(Md**2 - 1.0, 0.0)) * dlnV(Md)
    nu_d = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(Md))])
    nu = np.interp(M_a, Md, nu_d)

    cstar = pc / (rhot * at)
    return ExpansionTable(
        gamma_ref=float(gamma[0]),
        R=pc / (rhoc * Tc),
        pc=pc,
        Tc=Tc,
        M=M_a,
        p_pc=p_pc,
        T_Tc=T_a / Tc,
        rho_rhoc=rho_a / rhoc,
        a=a_a,
        gamma_local=gamma,
        nu=nu,
        area_ratio=eps_a,
        mode=f"cea_{mode}",
        meta={
            "cstar_m_s": float(cstar),
            "ox": ox,
            "fuel": fuel,
            "of": of,
            "pc_Pa": pc,
            "Tc_K": Tc,
            "eps_max": eps_max,
            "n": int(len(M_a)),
            "energy_residual_max": resid,  # |V_h/(M a) − 1|: CEA self-consistency
            "gamma_c": float(rhoc * ac**2 / pc),
            "gamma_t": float(gamma[0]),
            "gamma_e": float(gamma[-1]),
            "R_exit_over_R_chamber": float(p_abs[-1] / (rho_a[-1] * T_a[-1]) * rhoc * Tc / pc),
        },
    )


def cstar_cea(ox: str, fuel: str, of: float, pc: float) -> float:
    """CEA's own c* [m/s] (equilibrium) — used to cross-check the table."""
    from rocketcea.cea_obj import CEA_Obj

    return float(CEA_Obj(oxName=ox, fuelName=fuel).get_Cstar(pc / PSI, of)) * _FT


__all__ = ["MODES", "build_cea_table", "cstar_cea"]
