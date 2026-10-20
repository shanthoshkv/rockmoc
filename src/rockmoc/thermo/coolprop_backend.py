"""CoolProp backend — pure-fluid isentropic properties for cold-gas nozzles.

Optional: ``pip install rockmoc[coolprop]``.

For monatomic/diatomic cold-gas thrusters (N2, He, Ar, …), CoolProp provides
accurate real-gas properties. This backend builds a local effective-γ
PerfectGas surrogate at chamber conditions (ideal for preliminary MOC) and
exposes helpers for real isentropic stagnation relations along an expansion.
"""

from __future__ import annotations

from dataclasses import dataclass

from rockmoc.exceptions import ThermoError
from rockmoc.thermo.perfect_gas import PerfectGas


def _require_coolprop():
    try:
        import CoolProp.CoolProp as CP  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "CoolProp is not installed. Install with: pip install rockmoc[coolprop]"
        ) from exc
    return CP


@dataclass(slots=True)
class CoolPropGas:
    """Effective perfect-gas model from CoolProp fluid at (T0, P0)."""

    fluid: str = "Nitrogen"
    T0: float = 300.0  # K stagnation
    P0: float = 1.0e6  # Pa stagnation
    _gas: PerfectGas | None = None
    _summary: dict | None = None

    def solve(self) -> PerfectGas:
        CP = _require_coolprop()
        try:
            # cp, cv at stagnation state
            cp = float(CP.PropsSI("Cpmass", "T", self.T0, "P", self.P0, self.fluid))
            cv = float(CP.PropsSI("Cvmass", "T", self.T0, "P", self.P0, self.fluid))
            # Gas constant from EOS: R = cp - cv for ideal; for real use
            # specific gas constant approx from molar mass
            M = float(CP.PropsSI("M", "T", self.T0, "P", self.P0, self.fluid))  # kg/mol
            R_univ = 8.314462618  # J/(mol·K)
            R = R_univ / M  # J/(kg·K)
            gamma = cp / cv
            if gamma <= 1.0:
                raise ThermoError(f"CoolProp γ={gamma} invalid for {self.fluid}")
            # Prefer R from cp-cv if close to ideal
            R_alt = cp - cv
            if R_alt > 0 and abs(R_alt - R) / R < 0.2:
                R = R_alt
            pg = PerfectGas(gamma=gamma, R=R)
        except ThermoError:
            raise
        except Exception as exc:
            raise ThermoError(f"CoolProp evaluation failed: {exc}") from exc

        self._gas = pg
        self._summary = {
            "fluid": self.fluid,
            "T0_K": self.T0,
            "P0_Pa": self.P0,
            "gamma_eff": pg.gamma,
            "R": pg.R,
        }
        return pg

    @property
    def gas(self) -> PerfectGas:
        if self._gas is None:
            return self.solve()
        return self._gas

    @property
    def gamma(self) -> float:
        return self.gas.gamma

    @property
    def R(self) -> float:
        return self.gas.R

    def summary(self) -> dict:
        if self._summary is None:
            self.solve()
        assert self._summary is not None
        return dict(self._summary)

    def prandtl_meyer(self, M):
        return self.gas.prandtl_meyer(M)

    def prandtl_meyer_max(self):
        return self.gas.prandtl_meyer_max()

    def gamma_at_mach(self, M):
        return self.gas.gamma_at_mach(M)

    def mach_from_prandtl_meyer(self, nu):
        return self.gas.mach_from_prandtl_meyer(nu)

    def mach_angle(self, M):
        return self.gas.mach_angle(M)

    def area_ratio(self, M):
        return self.gas.area_ratio(M)

    def mach_from_area_ratio(self, area_ratio, *, branch="supersonic"):
        return self.gas.mach_from_area_ratio(area_ratio, branch=branch)

    def pressure_ratio(self, M):
        return self.gas.pressure_ratio(M)

    def temperature_ratio(self, M):
        return self.gas.temperature_ratio(M)
