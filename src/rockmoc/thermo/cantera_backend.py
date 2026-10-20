"""Cantera backend — equilibrium chamber → frozen γ for MOC design.

Optional: ``pip install rockmoc[cantera]``.

Computes HP equilibrium at chamber (P_c, H from reactants), freezes γ and
mean molecular weight at that state (or after isentropic expansion to sonic
as a simple throat estimate), and returns a :class:`PerfectGas` surrogate.
"""

from __future__ import annotations

from dataclasses import dataclass

from rockmoc.exceptions import ThermoError
from rockmoc.thermo.perfect_gas import PerfectGas


def _require_cantera():
    try:
        import cantera as ct  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "Cantera is not installed. Install with: pip install rockmoc[cantera]"
        ) from exc
    return ct


@dataclass(slots=True)
class CanteraGas:
    """Frozen gas model from a Cantera HP equilibrium chamber solve."""

    mechanism: str = "gri30.yaml"
    T_reactants: float = 300.0  # K
    P_chamber: float = 2.0e6  # Pa
    composition: str = "CH4:1, O2:2"
    _gas: PerfectGas | None = None
    _Tc: float | None = None
    _summary: dict | None = None

    def solve(self) -> PerfectGas:
        ct = _require_cantera()
        try:
            gas = ct.Solution(self.mechanism)
            gas.TPX = self.T_reactants, self.P_chamber, self.composition
            # Adiabatic flame / HP equilibrium at constant P
            gas.equilibrate("HP")
            Tc = float(gas.T)
            # Frozen γ = cp/cv at equilibrium composition (held fixed)
            gamma = float(gas.cp_mass / gas.cv_mass)
            # Mean molecular weight [kg/kmol] in Cantera is kg/kmol
            mw = float(gas.mean_molecular_weight)
            if gamma <= 1.0 or mw <= 0.0:
                raise ThermoError(f"Invalid Cantera state: γ={gamma}, MW={mw}")
            pg = PerfectGas.from_molecular_weight(gamma, mw)
        except ThermoError:
            raise
        except Exception as exc:
            raise ThermoError(f"Cantera evaluation failed: {exc}") from exc

        self._gas = pg
        self._Tc = Tc
        self._summary = {
            "mechanism": self.mechanism,
            "composition": self.composition,
            "P_chamber_Pa": self.P_chamber,
            "T_reactants_K": self.T_reactants,
            "Tc_K": Tc,
            "gamma_frozen": pg.gamma,
            "mw": mw,
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

    @property
    def Tc(self) -> float:
        if self._Tc is None:
            self.solve()
        assert self._Tc is not None
        return self._Tc

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
