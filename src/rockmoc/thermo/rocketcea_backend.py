"""RocketCEA backend — frozen / equilibrium chamber properties for design γ.

Optional dependency: ``pip install rockmoc[cea]`` (requires RocketCEA + CEA).

Industry use
------------
NASA CEA via RocketCEA is the de facto open path for rocket thermochemistry.
Typical design freezes composition at the throat and uses a single effective
γ for the nozzle MOC (frozen performance). This backend exposes that γ, MW,
T_c, and related scalars; full shifting-equilibrium tables are future work.
"""

from __future__ import annotations

from dataclasses import dataclass

from rockmoc.exceptions import ThermoError
from rockmoc.thermo.perfect_gas import PerfectGas


def _require_rocketcea():
    try:
        from rocketcea.cea_obj import CEA_Obj  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "RocketCEA is not installed. Install with: pip install rockmoc[cea]"
        ) from exc
    return CEA_Obj


@dataclass(slots=True)
class RocketCEAGas:
    """Frozen perfect-gas surrogate built from a RocketCEA chamber solution.

    Wraps :class:`PerfectGas` using γ and MW evaluated at the freeze point
    (default: throat), which is the standard engineering approach for MOC
    nozzle contouring with real propellants.
    """

    ox: str
    fuel: str
    of_ratio: float
    Pc_psia: float
    freeze_at: str = "throat"  # "throat" | "inject" | "chamber"
    eps: float = 40.0  # expansion ratio used for CEA call context
    _gas: PerfectGas | None = None
    _Tc: float | None = None
    _cstar: float | None = None
    _mw: float | None = None
    _gamma_throat: float | None = None
    _summary: dict | None = None

    def solve(self) -> PerfectGas:
        """Run CEA and return a PerfectGas with frozen throat γ and R."""
        CEA_Obj = _require_rocketcea()
        cea = CEA_Obj(oxName=self.ox, fuelName=self.fuel)
        # RocketCEA pressures in psia
        Pc = float(self.Pc_psia)
        ofr = float(self.of_ratio)
        eps = float(self.eps)

        try:
            # Molecular weight and gamma at chamber / throat
            mw, gamma = cea.get_Chamber_MolWt_gamma(Pc, ofr, eps)
            # Throat-focused values when available
            try:
                mw_t, gamma_t = cea.get_Throat_MolWt_gamma(Pc, ofr, eps)
            except Exception:
                mw_t, gamma_t = mw, gamma

            if self.freeze_at in ("throat", "t"):
                g_use, mw_use = float(gamma_t), float(mw_t)
            else:
                g_use, mw_use = float(gamma), float(mw)

            Tc = float(cea.get_Tcomb(Pc, ofr))
            try:
                cstar = float(cea.get_Cstar(Pc, ofr))
            except Exception:
                cstar = None

            isp_vac = None
            try:
                isp_vac = float(cea.estimate_Ambient_Isp(Pc, ofr, eps, Pamb=0.0)[0])
            except Exception:
                pass
        except Exception as exc:
            raise ThermoError(f"RocketCEA evaluation failed: {exc}") from exc

        gas = PerfectGas.from_molecular_weight(g_use, mw_use)
        self._gas = gas
        self._Tc = Tc
        self._cstar = cstar
        self._mw = mw_use
        self._gamma_throat = g_use
        self._summary = {
            "ox": self.ox,
            "fuel": self.fuel,
            "OF": ofr,
            "Pc_psia": Pc,
            "freeze_at": self.freeze_at,
            "gamma": g_use,
            "mw": mw_use,
            "Tc_K": Tc,
            "cstar_ft_s": cstar,
            "Isp_vac_s_est": isp_vac,
            "eps_context": eps,
        }
        return gas

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

    # Delegate MOC-facing methods to PerfectGas
    def prandtl_meyer(self, M):
        return self.gas.prandtl_meyer(M)

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
