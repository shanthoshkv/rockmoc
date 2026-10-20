"""High-fidelity axisymmetric MOC tests."""

from __future__ import annotations

import numpy as np

from rockmoc import DesignSpec, PerfectGas, design_ideal, design_mln, design_tic
from rockmoc.moc.axisym_hf import design_axisym_mln_hf, mesh_convergence_axisym


def test_mln_hf_basic():
    gas = PerfectGas(1.4)
    r = design_axisym_mln_hf(
        gas,
        throat_radius=0.05,
        exit_mach=2.5,
        n_chars=32,
        delta=1.0,
        match_exit_area=False,
        store_net=True,
    )
    assert r.kind == "mln"
    assert r.length > 0
    assert r.contour.r[-1] > 0.05
    assert 0.5 < r.length / 0.05 < 100
    assert np.all(np.isfinite(r.contour.x))
    s = r.summary()
    assert s["area_ratio_1d"] > 1.0
    assert s["area_ratio_geometric"] > 1.0


def test_public_design_mln_uses_hf():
    gas = PerfectGas(1.25)
    spec = DesignSpec(gas=gas, throat_radius=0.04, exit_mach=3.0, n_characteristics=30)
    r = design_mln(spec)
    assert r.kind == "mln"
    assert r.meta.get("designer") == "axisym_hf"
    assert r.metrics.length > 0


def test_ideal_and_tic():
    gas = PerfectGas(1.4)
    base = DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.4, n_characteristics=28)
    ideal = design_ideal(base)
    tic = design_tic(
        DesignSpec(
            gas=gas,
            throat_radius=0.05,
            exit_mach=2.4,
            n_characteristics=28,
            truncate_fraction=0.7,
        )
    )
    assert tic.metrics.length < ideal.metrics.length
    assert ideal.metrics.length > 0


def test_planar_delta0():
    gas = PerfectGas(1.4)
    r = design_axisym_mln_hf(
        gas,
        throat_radius=0.05,
        exit_mach=2.0,
        n_chars=24,
        delta=0.0,
        match_exit_area=False,
    )
    assert r.length > 0
    # θ_max = νe/2 for planar
    nu = float(gas.prandtl_meyer(2.0))
    assert abs(r.theta_max_rad - 0.5 * nu) < 1e-9


def test_mesh_convergence_runs():
    gas = PerfectGas(1.4)
    rows = mesh_convergence_axisym(
        gas, throat_radius=0.05, exit_mach=2.2, n_list=[16, 24], delta=1.0
    )
    assert len(rows) == 2
    assert rows[0]["length_over_Rt"] > 0


def test_match_exit_area_option():
    gas = PerfectGas(1.4)
    Me = 2.5
    eps = float(gas.area_ratio(Me))
    r = design_axisym_mln_hf(
        gas,
        throat_radius=0.05,
        exit_mach=Me,
        n_chars=24,
        match_exit_area=True,
        store_net=False,
    )
    assert abs(r.area_ratio_geometric - eps) / eps < 1e-6
