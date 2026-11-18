"""Real-gas (RocketCEA equilibrium/frozen) expansion table: accuracy and wiring."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("rocketcea")

from rocketcea.cea_obj import CEA_Obj  # noqa: E402

from rockmoc.designers import DesignSpec  # noqa: E402
from rockmoc.exceptions import DesignError  # noqa: E402
from rockmoc.spec_io import build_design, performance_from_spec, resolve_gas  # noqa: E402
from rockmoc.thermo import VariableGammaGas  # noqa: E402
from rockmoc.thermo.cea_table import PSI, build_cea_table, cstar_cea  # noqa: E402

PC, OF = 7e6, 2.4
G0 = 9.80665


@pytest.mark.parametrize("mode", ["equilibrium", "frozen"])
def test_table_reproduces_cea_isp_mach_cstar(mode):
    tb = build_cea_table("LOX", "RP1", OF, PC, mode=mode)
    g = VariableGammaGas(tb)
    cea = CEA_Obj(oxName="LOX", fuelName="RP1")
    fz = int(mode == "frozen")
    for eps in (5.0, 25.0, 100.0):
        Me = float(g.mach_from_area_ratio(eps))
        ref = cea.get_MachNumber(PC / PSI, OF, eps, frozen=fz, frozenAtThroat=fz)
        assert Me == pytest.approx(ref, rel=1e-3)
        cf = g.thrust_coefficient(float(g.pe_pc_from_mach(Me)), 0.0, exit_mach=Me, area_ratio=eps)
        isp = cf * g.characteristic_velocity(tb.Tc) / G0
        ivac = cea.get_IvacCstrTc(PC / PSI, OF, eps, frozen=fz, frozenAtThroat=fz)[0]
        assert isp == pytest.approx(ivac, rel=1e-3)
    assert tb.meta["cstar_m_s"] == pytest.approx(cstar_cea("LOX", "RP1", OF, PC), rel=1e-4)
    assert tb.meta["energy_residual_max"] < 1e-3  # V from h vs M·a: CEA self-consistent
    assert tb.gamma_local[-1] > tb.gamma_local[0]  # γ rises as the gas cools


def test_independent_cantera_expansion():
    ct = pytest.importorskip("cantera")
    tb = build_cea_table("LOX", "LH2", 6.0, PC, mode="equilibrium", eps_max=200)
    g = VariableGammaGas(tb)
    gas = ct.Solution("h2o2.yaml")
    gas.TPY = tb.Tc, PC, "O2:6, H2:1"
    gas.equilibrate("TP")
    hc, sc = gas.enthalpy_mass, gas.entropy_mass
    for eps in (5.0, 20.0, 80.0):
        Me = float(g.mach_from_area_ratio(eps))
        gas.SP = sc, float(g.pe_pc_from_mach(Me)) * PC
        gas.equilibrate("SP")
        Te = tb.Tc / float(g.temperature_ratio(Me))
        Ve = float(np.sqrt(2.0 * (hc - gas.enthalpy_mass)))
        assert Te == pytest.approx(gas.T, rel=0.01)
        assert Ve == pytest.approx(Me * float(np.interp(Me, tb.M, tb.a)), rel=0.005)


SPEC = {
    "design": {"type": "mln", "area_ratio": 12.0, "n_characteristics": 20},
    "propellants": {"ox": "LOX", "fuel": "RP1", "of": OF, "pc": PC},
    "performance": {"pc": PC},
}


def test_design_uses_real_gas_end_to_end():
    gas, Tc = resolve_gas(SPEC)  # default model with [propellants] = equilibrium
    assert isinstance(gas, VariableGammaGas) and Tc == pytest.approx(gas.table.Tc)
    res = build_design(SPEC)
    assert res.metrics.area_ratio == pytest.approx(12.0, rel=5e-3)
    perf = performance_from_spec(res, SPEC)
    ivac = CEA_Obj(oxName="LOX", fuelName="RP1").get_IvacCstrTc(PC / PSI, OF, 12.0)[0]
    # ideal 1-D Cf·c*/g0 must equal CEA; the loss stack (divergence, exit KE) only reduces it
    assert perf["Cf_1d_vac"] * perf["cstar_m_s"] / G0 == pytest.approx(ivac, rel=1e-3)
    assert 0.97 * ivac < perf["Isp_vac_s"] <= ivac * 1.001
    frozen = build_design({**SPEC, "gas": {"model": "frozen"}})
    assert frozen.metrics.area_ratio == pytest.approx(12.0, rel=5e-3)
    assert performance_from_spec(frozen, {**SPEC, "gas": {"model": "frozen"}})["Isp_vac_s"] < perf[
        "Isp_vac_s"
    ]  # frozen flow loses recombination energy


def test_pressure_ratio_inverts_real_gas_and_range_guard():
    gas, _ = resolve_gas(SPEC)
    ds = DesignSpec(gas=gas, pressure_ratio=0.01)
    assert float(gas.pe_pc_from_mach(ds.resolved_exit_mach())) == pytest.approx(0.01, rel=1e-6)
    with pytest.raises(DesignError):
        DesignSpec(gas=gas, area_ratio=1e5).resolved_exit_mach()
    with pytest.raises(DesignError):
        DesignSpec(gas=gas, pressure_ratio=1e-9).resolved_exit_mach()
    with pytest.raises(DesignError):  # real-gas model without propellants
        resolve_gas({"gas": {"model": "frozen"}})
