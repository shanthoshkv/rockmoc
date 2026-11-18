"""Isentropic flow states from stagnation, for constant-γ and real-gas (table) models."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def edge_state(
    gas: object, M: ArrayLike, Tc: float, pc: float
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """``(T, p, rho, a)`` at Mach ``M`` on the isentrope from stagnation ``(Tc, pc)``.

    Perfect gas: closed form with ``gas.gamma``. Real-gas model (has ``.table``): the
    table's own T, p, ρ, a (γ(T), composition vary) for M ≥ 1, scaled to the requested
    ``Tc``/``pc``; the subsonic part (M < 1) uses the throat-γ closed form.
    """
    M = np.asarray(M, dtype=np.float64)
    g, R = float(gas.gamma), float(gas.R)  # type: ignore[attr-defined]
    m2 = 0.5 * (g - 1.0) * M**2
    T = Tc / (1.0 + m2)
    p = pc / (1.0 + m2) ** (g / (g - 1.0))
    rho = p / (R * T)
    a = np.sqrt(g * R * T)
    tab = getattr(gas, "table", None)
    if tab is None:
        return T, p, rho, a
    sup = M >= 1.0
    Ms = np.clip(M, tab.M[0], tab.M[-1])
    T_t = Tc / np.asarray(gas.temperature_ratio(Ms))  # type: ignore[attr-defined]
    p_t = pc / np.asarray(gas.pressure_ratio(Ms))  # type: ignore[attr-defined]
    rho_t = np.interp(Ms, tab.M, tab.rho_rhoc) * pc / (tab.R * Tc)
    a_t = np.interp(Ms, tab.M, tab.a) * np.sqrt(Tc / tab.Tc)
    return (np.where(sup, T_t, T), np.where(sup, p_t, p),
            np.where(sup, rho_t, rho), np.where(sup, a_t, a))
