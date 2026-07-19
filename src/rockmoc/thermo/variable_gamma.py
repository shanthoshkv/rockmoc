"""Variable-γ gas adapter: ExpansionTable → PerfectGas-like API for MOC.

Uses local γ(M) and integrated ν(M) from the table for Prandtl–Meyer and
area–Mach, while isentropic p,T ratios use local γ at the evaluation point.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from rockmoc.exceptions import ThermoError
from rockmoc.thermo.expansion_table import ExpansionTable
from rockmoc.thermo.perfect_gas import PerfectGas
from rockmoc.thermo.solve import bracketed_newton


@dataclass
class VariableGammaGas:
    """MOC-facing gas model backed by an :class:`ExpansionTable`."""

    table: ExpansionTable

    @property
    def gamma(self) -> float:
        """Reference γ (throat / design) for APIs expecting a scalar."""
        return float(self.table.gamma_ref)

    @property
    def R(self) -> float:
        return float(self.table.R)

    def gamma_at_mach(self, M: float) -> float:
        M = float(np.clip(M, self.table.M[0], self.table.M[-1]))
        return float(np.interp(M, self.table.M, self.table.gamma_local))

    def prandtl_meyer(self, M: ArrayLike) -> NDArray | float:
        arr = np.asarray(M, dtype=np.float64)
        out = np.interp(arr, self.table.M, self.table.nu)
        return float(out) if out.ndim == 0 else out

    def mach_from_prandtl_meyer(self, nu: ArrayLike) -> NDArray | float:
        arr = np.asarray(nu, dtype=np.float64)
        out = np.interp(arr, self.table.nu, self.table.M)
        return float(out) if np.ndim(out) == 0 else out

    def mach_angle(self, M: ArrayLike) -> NDArray | float:
        arr = np.asarray(M, dtype=np.float64)
        if np.any(arr < 1.0):
            raise ThermoError("Mach angle requires M ≥ 1")
        out = np.arcsin(1.0 / arr)
        return float(out) if out.ndim == 0 else out

    def area_ratio(self, M: ArrayLike) -> NDArray | float:
        arr = np.asarray(M, dtype=np.float64)
        out = np.interp(arr, self.table.M, self.table.area_ratio)
        return float(out) if out.ndim == 0 else out

    def mach_from_area_ratio(
        self, area_ratio: ArrayLike, *, branch: str = "supersonic"
    ) -> NDArray | float:
        if branch != "supersonic":
            # fall back to perfect gas for subsonic
            return PerfectGas(self.gamma, self.R).mach_from_area_ratio(
                area_ratio, branch=branch
            )
        ar = np.asarray(area_ratio, dtype=np.float64)
        out = np.interp(ar, self.table.area_ratio, self.table.M)
        return float(out) if np.ndim(out) == 0 else out

    def pressure_ratio(self, M: ArrayLike) -> NDArray | float:
        """p0/p."""
        arr = np.asarray(M, dtype=np.float64)
        p_pc = np.interp(arr, self.table.M, self.table.p_pc)
        out = 1.0 / np.maximum(p_pc, 1e-30)
        return float(out) if out.ndim == 0 else out

    def temperature_ratio(self, M: ArrayLike) -> NDArray | float:
        """T0/T."""
        arr = np.asarray(M, dtype=np.float64)
        T_Tc = np.interp(arr, self.table.M, self.table.T_Tc)
        out = 1.0 / np.maximum(T_Tc, 1e-30)
        return float(out) if out.ndim == 0 else out

    def pe_pc_from_mach(self, M: ArrayLike) -> NDArray | float:
        arr = np.asarray(M, dtype=np.float64)
        out = np.interp(arr, self.table.M, self.table.p_pc)
        return float(out) if out.ndim == 0 else out

    def prandtl_meyer_max(self) -> float:
        return float(self.table.nu[-1])

    def characteristic_velocity(self, T0: float) -> float:
        g = self.gamma
        return float(
            np.sqrt(self.R * T0 / g)
            * ((g + 1.0) / 2.0) ** (0.5 * (g + 1.0) / (g - 1.0))
        )

    def thrust_coefficient(
        self,
        pe_pc: float,
        pa_pc: float,
        *,
        area_ratio: float | None = None,
        exit_mach: float | None = None,
    ) -> float:
        g = self.gamma
        if exit_mach is None:
            if area_ratio is None:
                raise ThermoError("Provide area_ratio or exit_mach")
            Me = float(self.mach_from_area_ratio(area_ratio, branch="supersonic"))
        else:
            Me = float(exit_mach)
            if area_ratio is None:
                area_ratio = float(self.area_ratio(Me))
        term = (2.0 * g**2 / (g - 1.0)) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0))
        pe_term = 1.0 - pe_pc ** ((g - 1.0) / g)
        return float(np.sqrt(term * pe_term) + area_ratio * (pe_pc - pa_pc))

    def as_perfect_gas(self) -> PerfectGas:
        return PerfectGas(gamma=self.gamma, R=self.R)

    def gamma_along_length(
        self, x: ArrayLike, x_throat: float = 0.0, x_exit: float | None = None
    ) -> NDArray | float:
        """Map axial station to γ via local M(x) schedule (linear M in x).

        Used for post-processing γ(x). MOC itself uses γ(M) via the table.
        """
        arr = np.asarray(x, dtype=np.float64)
        if x_exit is None:
            x_exit = float(np.max(arr)) if arr.size else 1.0
        span = max(float(x_exit) - float(x_throat), 1e-12)
        # Map x → M linearly between table ends for a smooth γ(x)
        M = self.table.M[0] + (self.table.M[-1] - self.table.M[0]) * (
            (arr - x_throat) / span
        )
        M = np.clip(M, self.table.M[0], self.table.M[-1])
        out = np.interp(M, self.table.M, self.table.gamma_local)
        return float(out) if out.ndim == 0 else out

    def local_gamma_from_nu(self, nu: float) -> float:
        """γ at the Mach corresponding to Prandtl–Meyer angle ν."""
        M = float(self.mach_from_prandtl_meyer(nu))
        return self.gamma_at_mach(M)
