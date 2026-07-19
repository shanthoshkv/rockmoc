"""Perfect-gas isentropic relations and Prandtl–Meyer functions.

All angles returned by Prandtl–Meyer helpers are in **radians**.
All quantities use SI units unless noted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from rockmoc.exceptions import ThermoError
from rockmoc.thermo.solve import bracketed_newton

FloatArray = NDArray[np.floating]


def _as_float_array(x: ArrayLike) -> FloatArray:
    return np.asarray(x, dtype=np.float64)


def _validate_gamma(gamma: float) -> float:
    g = float(gamma)
    if not (1.05 < g <= 5.0 / 3.0 + 1e-12):
        # Allow slightly above 5/3 for numerical convenience; hard reject nonsense.
        if g <= 1.0 or g > 2.0:
            raise ThermoError(f"gamma must be > 1 and typically ≤ 5/3; got {g}")
    return g


def _validate_mach(M: ArrayLike, *, allow_subsonic: bool = True) -> FloatArray:
    arr = _as_float_array(M)
    if np.any(arr <= 0.0):
        raise ThermoError("Mach number must be positive")
    if not allow_subsonic and np.any(arr < 1.0):
        raise ThermoError("Mach number must be ≥ 1 for this relation")
    return arr


@dataclass(frozen=True, slots=True)
class PerfectGas:
    """Constant-γ perfect gas.

    Parameters
    ----------
    gamma:
        Ratio of specific heats (dimensionless).
    R:
        Specific gas constant [J/(kg·K)]. If only molecular weight is known,
        use :meth:`from_molecular_weight`.
    """

    gamma: float = 1.4
    R: float = 287.052874  # dry air default, J/(kg·K)

    def __post_init__(self) -> None:
        object.__setattr__(self, "gamma", _validate_gamma(self.gamma))
        if self.R <= 0.0:
            raise ThermoError(f"R must be positive; got {self.R}")

    @classmethod
    def from_molecular_weight(
        cls,
        gamma: float,
        molecular_weight: float,
        *,
        R_univ: float = 8314.462618,  # J/(kmol·K)
    ) -> PerfectGas:
        """Build gas from γ and molecular weight [kg/kmol] (or g/mol)."""
        if molecular_weight <= 0.0:
            raise ThermoError("molecular_weight must be positive")
        return cls(gamma=gamma, R=R_univ / molecular_weight)

    @property
    def cp(self) -> float:
        """Specific heat at constant pressure [J/(kg·K)]."""
        return self.R * self.gamma / (self.gamma - 1.0)

    @property
    def cv(self) -> float:
        """Specific heat at constant volume [J/(kg·K)]."""
        return self.R / (self.gamma - 1.0)

    # ------------------------------------------------------------------
    # Isentropic relations
    # ------------------------------------------------------------------

    def temperature_ratio(self, M: ArrayLike) -> FloatArray:
        """T0/T for given Mach number(s)."""
        M = _validate_mach(M)
        g = self.gamma
        return 1.0 + 0.5 * (g - 1.0) * M**2

    def pressure_ratio(self, M: ArrayLike) -> FloatArray:
        """p0/p for given Mach number(s)."""
        g = self.gamma
        return self.temperature_ratio(M) ** (g / (g - 1.0))

    def density_ratio(self, M: ArrayLike) -> FloatArray:
        """ρ0/ρ for given Mach number(s)."""
        g = self.gamma
        return self.temperature_ratio(M) ** (1.0 / (g - 1.0))

    def area_ratio(self, M: ArrayLike) -> FloatArray:
        """A/A* for given Mach number(s) (isentropic)."""
        M = _validate_mach(M)
        g = self.gamma
        expo = 0.5 * (g + 1.0) / (g - 1.0)
        term = (2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M**2)
        return (1.0 / M) * term**expo

    def mach_from_area_ratio(
        self,
        area_ratio: ArrayLike,
        *,
        branch: str = "supersonic",
    ) -> FloatArray:
        """Invert A/A* → M on the subsonic or supersonic branch."""
        ar = _as_float_array(area_ratio)
        if np.any(ar < 1.0 - 1e-14):
            raise ThermoError("area ratio A/A* must be ≥ 1")
        if branch not in {"supersonic", "subsonic"}:
            raise ThermoError("branch must be 'supersonic' or 'subsonic'")

        g = self.gamma
        out = np.empty_like(ar, dtype=np.float64)

        def solve_one(eps: float) -> float:
            if abs(eps - 1.0) < 1e-14:
                return 1.0

            def f(M: float) -> float:
                return float(self.area_ratio(M) - eps)

            if branch == "supersonic":
                # A/A*(M) is minimized at M=1 and increases for M>1.
                # f = A/A* - eps ≤ 0 near M=1 (eps≥1) and becomes >0 for large M.
                lo, hi = 1.0 + 1e-12, 2.0
                while f(hi) < 0.0 and hi < 1.0e4:
                    hi *= 1.5
                if f(hi) < 0.0:
                    raise ThermoError(f"Could not bracket supersonic M for A/A*={eps}")
                x0 = min(max(1.5, lo), hi)
                return bracketed_newton(f, None, x0, x_min=lo, x_max=hi)
            else:
                # Subsonic: A/A* decreases from ∞ as M→0 toward 1 at M=1.
                lo, hi = 1e-8, 1.0 - 1e-12
                # Ensure f(lo) and f(hi) have opposite signs when possible.
                if f(lo) < 0.0:
                    # eps large: need smaller M
                    while f(lo) < 0.0 and lo > 1e-12:
                        lo *= 0.5
                return bracketed_newton(f, None, 0.5 * (lo + hi), x_min=lo, x_max=hi)

        it = np.nditer(
            [ar, out],
            flags=["multi_index", "refs_ok"],
            op_flags=[["readonly"], ["writeonly"]],
        )
        for a, o in it:
            o[...] = solve_one(float(a))
        return out if out.shape else float(out)

    def prandtl_meyer(self, M: ArrayLike) -> FloatArray:
        """Prandtl–Meyer function ν(M) [radians]."""
        M = _validate_mach(M, allow_subsonic=False)
        g = self.gamma
        gp = np.sqrt((g + 1.0) / (g - 1.0))
        # ν = gp * arctan(sqrt((M^2-1)/gp^2 wait)) standard form:
        # ν = sqrt((γ+1)/(γ-1)) * arctan(sqrt(((γ-1)/(γ+1))*(M^2-1))) - arctan(sqrt(M^2-1))
        root = np.sqrt(np.maximum(M**2 - 1.0, 0.0))
        return gp * np.arctan(root / gp) - np.arctan(root)

    def mach_from_prandtl_meyer(self, nu: ArrayLike) -> FloatArray:
        """Invert ν(M) → M (supersonic). ν in radians."""
        nu_arr = _as_float_array(nu)
        if np.any(nu_arr < -1e-15):
            raise ThermoError("Prandtl–Meyer angle must be ≥ 0")

        g = self.gamma
        nu_max = float(self.prandtl_meyer_max())
        if np.any(nu_arr > nu_max + 1e-12):
            raise ThermoError(
                f"ν exceeds maximum for γ={g}: ν_max={nu_max:.6f} rad "
                f"({np.degrees(nu_max):.4f} deg)"
            )

        out = np.empty_like(nu_arr, dtype=np.float64)

        def solve_one(target: float) -> float:
            if target <= 1e-15:
                return 1.0

            def f(M: float) -> float:
                return float(self.prandtl_meyer(M) - target)

            # dν/dM = sqrt(M^2-1)/M * 1/(1 + 0.5(γ-1)M^2) * ... analytical optional
            lo, hi = 1.0 + 1e-12, 50.0
            # Ensure upper bound is high enough (ν is monotonic).
            while f(hi) < 0.0 and hi < 1e4:
                hi *= 1.5
            return bracketed_newton(f, None, 2.0, x_min=lo, x_max=hi)

        for idx, val in np.ndenumerate(nu_arr):
            out[idx] = solve_one(float(val))
        return out if out.shape else float(out)

    def prandtl_meyer_max(self) -> float:
        """Maximum Prandtl–Meyer angle [radians] as M→∞."""
        g = self.gamma
        gp = np.sqrt((g + 1.0) / (g - 1.0))
        return float(0.5 * np.pi * (gp - 1.0))

    def mach_angle(self, M: ArrayLike) -> FloatArray:
        """Mach angle μ = arcsin(1/M) [radians]."""
        M = _validate_mach(M, allow_subsonic=False)
        return np.arcsin(1.0 / M)

    def sound_speed(self, T: ArrayLike) -> FloatArray:
        """Local speed of sound a = sqrt(γ R T) [m/s]."""
        T_arr = _as_float_array(T)
        if np.any(T_arr <= 0.0):
            raise ThermoError("Temperature must be positive")
        return np.sqrt(self.gamma * self.R * T_arr)

    def velocity(self, M: ArrayLike, T: ArrayLike) -> FloatArray:
        """Flow speed V = M * a(T) [m/s]."""
        return _validate_mach(M) * self.sound_speed(T)

    def characteristic_velocity(self, T0: float) -> float:
        """Ideal characteristic velocity c* [m/s] for chamber temperature T0."""
        if T0 <= 0.0:
            raise ThermoError("T0 must be positive")
        g = self.gamma
        return float(
            np.sqrt(self.R * T0 / g) * ((g + 1.0) / 2.0) ** (0.5 * (g + 1.0) / (g - 1.0))
        )

    def thrust_coefficient(
        self,
        pe_pc: float,
        pa_pc: float,
        *,
        area_ratio: float | None = None,
        exit_mach: float | None = None,
    ) -> float:
        """Ideal 1D thrust coefficient Cf.

        Provide either ``area_ratio`` or ``exit_mach``. ``pe_pc`` is exit static
        pressure over chamber pressure; ``pa_pc`` is ambient over chamber.
        """
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
        cf = np.sqrt(term * pe_term) + area_ratio * (pe_pc - pa_pc)
        return float(cf)

    def pe_pc_from_mach(self, M: ArrayLike) -> FloatArray:
        """Static exit pressure ratio pe/p0 (chamber) for isentropic expansion."""
        return 1.0 / self.pressure_ratio(M)


# Module-level convenience functions (default air γ=1.4) ---------------------------------

_DEFAULT = PerfectGas()


def area_ratio(M: ArrayLike, gamma: float = 1.4) -> FloatArray:
    return PerfectGas(gamma=gamma).area_ratio(M)


def mach_from_area_ratio(
    area_ratio_val: ArrayLike,
    gamma: float = 1.4,
    *,
    branch: str = "supersonic",
) -> FloatArray:
    return PerfectGas(gamma=gamma).mach_from_area_ratio(area_ratio_val, branch=branch)


def prandtl_meyer(M: ArrayLike, gamma: float = 1.4) -> FloatArray:
    return PerfectGas(gamma=gamma).prandtl_meyer(M)


def mach_from_prandtl_meyer(nu: ArrayLike, gamma: float = 1.4) -> FloatArray:
    return PerfectGas(gamma=gamma).mach_from_prandtl_meyer(nu)


def mach_angle(M: ArrayLike) -> FloatArray:
    return _DEFAULT.mach_angle(M)
