"""Variable-γ gas adapter: ExpansionTable → PerfectGas-like API for MOC.

Uses local γ(M) and integrated ν(M) from the table for Prandtl–Meyer and
area–Mach, while isentropic p,T ratios use local γ at the evaluation point.
"""

from __future__ import annotations

import math
from bisect import bisect_right
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import PchipInterpolator

from rockmoc.exceptions import ThermoError
from rockmoc.thermo.expansion_table import ExpansionTable
from rockmoc.thermo.perfect_gas import PerfectGas


@dataclass
class VariableGammaGas:
    """MOC-facing gas model backed by an :class:`ExpansionTable`."""

    table: ExpansionTable

    def __post_init__(self) -> None:
        # Monotone cubic (PCHIP) interpolants: 3rd-order accurate and shape-
        # preserving, vs. 1st-order np.interp (matters for ν(M) in MOC).
        t = self.table
        self._nu_M = PchipInterpolator(t.M, t.nu, extrapolate=False)
        self._M_nu = PchipInterpolator(t.nu, t.M, extrapolate=False)
        self._ar_M = PchipInterpolator(t.M, t.area_ratio, extrapolate=False)
        self._M_ar = PchipInterpolator(t.area_ratio, t.M, extrapolate=False)
        # Dense samples of the PCHIP curves for fast scalar lookups (MOC calls
        # these per point): linear on 4000 nodes, error ~1e-7 relative.
        Md = np.linspace(t.M[0], t.M[-1], 4000)
        nud = self._nu_M(Md)
        ard = self._ar_M(Md)
        gd = np.interp(Md, t.M, t.gamma_local)
        self._dense = {k: v.tolist() for k, v in (("M", Md), ("nu", nud), ("ar", ard), ("g", gd))}

    def _lookup(self, xkey: str, ykey: str, x: float) -> float:
        xs, ys = self._dense[xkey], self._dense[ykey]
        if x <= xs[0]:
            return ys[0]
        if x >= xs[-1]:
            return ys[-1]
        i = bisect_right(xs, x) - 1
        f = (x - xs[i]) / (xs[i + 1] - xs[i])
        return ys[i] + f * (ys[i + 1] - ys[i])

    @staticmethod
    def _eval(f: PchipInterpolator, x: ArrayLike, lo: float, hi: float) -> NDArray | float:
        out = f(np.clip(np.asarray(x, dtype=np.float64), lo, hi))
        return float(out) if np.ndim(out) == 0 else out

    @property
    def gamma(self) -> float:
        """Reference γ (throat / design) for APIs expecting a scalar."""
        return float(self.table.gamma_ref)

    @property
    def R(self) -> float:
        return float(self.table.R)

    def gamma_at_mach(self, M: float) -> float:
        if isinstance(M, float):
            return self._lookup("M", "g", M)
        M = float(np.clip(M, self.table.M[0], self.table.M[-1]))
        return float(np.interp(M, self.table.M, self.table.gamma_local))

    def prandtl_meyer(self, M: ArrayLike) -> NDArray | float:
        if isinstance(M, float):
            return self._lookup("M", "nu", M)
        return self._eval(self._nu_M, M, self.table.M[0], self.table.M[-1])

    def mach_from_prandtl_meyer(self, nu: ArrayLike) -> NDArray | float:
        if isinstance(nu, float):
            return self._lookup("nu", "M", nu)
        return self._eval(self._M_nu, nu, self.table.nu[0], self.table.nu[-1])

    def mach_angle(self, M: ArrayLike) -> NDArray | float:
        if isinstance(M, float) and M >= 1.0:
            return math.asin(1.0 / M)
        arr = np.asarray(M, dtype=np.float64)
        if np.any(arr < 1.0):
            raise ThermoError("Mach angle requires M ≥ 1")
        out = np.arcsin(1.0 / arr)
        return float(out) if out.ndim == 0 else out

    def area_ratio(self, M: ArrayLike) -> NDArray | float:
        if isinstance(M, float):
            return self._lookup("M", "ar", M)
        return self._eval(self._ar_M, M, self.table.M[0], self.table.M[-1])

    def mach_from_area_ratio(
        self, area_ratio: ArrayLike, *, branch: str = "supersonic"
    ) -> NDArray | float:
        if branch != "supersonic":
            # fall back to perfect gas for subsonic
            return PerfectGas(self.gamma, self.R).mach_from_area_ratio(
                area_ratio, branch=branch
            )
        return self._eval(self._M_ar, area_ratio, self.table.area_ratio[0], self.table.area_ratio[-1])

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
        cs = self.table.meta.get("cstar_m_s")
        if cs is not None:  # integrated isentrope: ratios depend on γ(M) only, c* ∝ √T0
            return float(cs) * math.sqrt(T0 / self.table.Tc)
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
        """Cf = ṁV_e/(p_c A*) + ε (p_e − p_a)/p_c with the table's own state.

        ρV² = γ p M² for any ideal gas, so the momentum term is
        (A_e/A*)(M_e) · (p_e/p_c) · γ_e M_e² — exact for variable γ.
        """
        if exit_mach is None:
            if area_ratio is None:
                raise ThermoError("Provide area_ratio or exit_mach")
            Me = float(self.mach_from_area_ratio(area_ratio, branch="supersonic"))
        else:
            Me = float(exit_mach)
            if area_ratio is None:
                area_ratio = float(self.area_ratio(Me))
        pe_tab = float(self.pe_pc_from_mach(Me))
        momentum = float(self.area_ratio(Me)) * pe_tab * self.gamma_at_mach(Me) * Me * Me
        return float(momentum + area_ratio * (pe_pc - pa_pc))

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
