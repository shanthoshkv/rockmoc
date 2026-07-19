"""Optional thermo backends — skip cleanly if deps missing."""

from __future__ import annotations

import pytest


def test_rocketcea_optional():
    pytest.importorskip("rocketcea")
    from rockmoc.thermo.rocketcea_backend import RocketCEAGas

    # LOX/LH2 style; Pc in psia for RocketCEA
    backend = RocketCEAGas(ox="LOX", fuel="LH2", of_ratio=5.0, Pc_psia=500.0, eps=20.0)
    gas = backend.solve()
    assert gas.gamma > 1.0
    assert gas.R > 0.0
    s = backend.summary()
    assert s["Tc_K"] > 1000.0


def test_cantera_optional():
    pytest.importorskip("cantera")
    from rockmoc.thermo.cantera_backend import CanteraGas

    backend = CanteraGas(
        mechanism="gri30.yaml",
        composition="CH4:1, O2:2",
        P_chamber=1.0e6,
        T_reactants=300.0,
    )
    gas = backend.solve()
    assert 1.05 < gas.gamma < 1.7
    assert backend.Tc > 1000.0


def test_coolprop_optional():
    pytest.importorskip("CoolProp")
    from rockmoc.thermo.coolprop_backend import CoolPropGas

    backend = CoolPropGas(fluid="Nitrogen", T0=300.0, P0=2.0e6)
    gas = backend.solve()
    assert 1.2 < gas.gamma < 1.5
    assert gas.R > 200.0
