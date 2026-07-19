"""Tests for production pipeline A–F (no CFD)."""

from __future__ import annotations

import math

import pytest

from rockmoc import (
    DesignSpec,
    PerfectGas,
    design_mln,
    design_production,
    design_rao_toc,
    field_integrated_cf,
    mesh_certificate,
    monte_carlo_mln,
    sauer_throat_ivl,
    sensitivity_one_at_a_time,
)
from rockmoc.thermo import VariableGammaGas, build_multi_gamma_table, build_perfect_gas_table


def test_production_pipeline_mln():
    gas = PerfectGas(1.25)
    spec = DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.8, n_characteristics=24)
    # n_chars_min lowered for test speed
    prod = design_production(
        spec,
        kind="mln",
        apply_bl=True,
        n_chars_min=24,
        certificate=True,
        store_net=True,
        altitude_points=[0.0, 20e3],
    )
    s = prod.summary()
    assert prod.design.kind == "mln"
    assert prod.field.Cf_vac > 0.5
    assert prod.bl is not None
    assert prod.bl.delta_star[-1] > 0
    assert prod.table is not None
    assert "M_axis" in prod.throat_ivl_meta
    assert len(prod.altitude_sweep) == 2
    assert "production" in s


def test_production_multi_gamma():
    gas = PerfectGas(1.22)
    spec = DesignSpec(gas=gas, throat_radius=0.04, exit_mach=2.5, n_characteristics=20)
    prod = design_production(
        spec,
        multi_gamma=(1.25, 1.22, 1.18),
        apply_bl=False,
        n_chars_min=20,
        certificate=False,
        store_net=False,
    )
    assert prod.table is not None
    assert prod.table.mode == "multi_gamma"
    assert prod.field.Cf_vac > 0.5


def test_field_cf():
    gas = PerfectGas(1.4)
    r = design_mln(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.4, n_characteristics=20),
        store_net=True,
    )
    f = field_integrated_cf(
        r.contour,
        gas,
        throat_radius=0.05,
        exit_mach_design=2.4,
        pc=5e6,
        Tc=3000.0,
        pa=0.0,
        net=r.net,
    )
    assert f.Cf_vac > 0.5
    assert f.mass_flow > 0 or f.method == "contour_exit_station"
    assert math.isfinite(f.thrust_N)


def test_variable_gamma_gas():
    tab = build_multi_gamma_table(
        gamma_c=1.25, gamma_t=1.22, gamma_e=1.18, R=350.0, M_max=5.0, n=60
    )
    vg = VariableGammaGas(tab)
    assert vg.gamma == pytest.approx(1.22)
    nu = float(vg.prandtl_meyer(2.5))
    M2 = float(vg.mach_from_prandtl_meyer(nu))
    assert abs(M2 - 2.5) < 0.05
    assert float(vg.area_ratio(2.0)) > 1.0


def test_rao_toc():
    gas = PerfectGas(1.2)
    r = design_rao_toc(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=20),
        length_fraction=0.75,
    )
    assert r.kind == "rao"
    assert r.meta.get("method") == "rao_toc_moc_ideal_length"
    assert r.metrics.length > 0
    # Should be shorter than ideal
    assert r.meta["length_m"] < r.meta["L_ideal_moc_m"] * 0.95


def test_sauer_hall_ivl():
    gas = PerfectGas(1.3)
    ivl = sauer_throat_ivl(gas, throat_radius=0.05, rc_over_rt=1.5, n_points=21)
    assert float(ivl.M.min()) > 1.0
    assert "hall" in str(ivl.meta.get("model", "sauer")).lower() or ivl.meta.get("A", 0) > 0


def test_uq_and_certificate():
    uq = monte_carlo_mln(
        throat_radius=0.05,
        exit_mach=2.5,
        gamma_mean=1.25,
        gamma_std=0.015,
        n_chars_choices=[16, 24],
        n_samples=6,
        seed=1,
    )
    assert uq.stats["L_over_Rt"]["std"] >= 0
    assert uq.meta["n_ok"] >= 3
    sens = sensitivity_one_at_a_time(
        throat_radius=0.05, exit_mach=2.5, gamma=1.25, n_chars=16, d_gamma=0.02
    )
    assert "dL_dgamma" in sens
    cert = mesh_certificate(
        throat_radius=0.05, exit_mach=2.4, gamma=1.4, n_coarse=16, n_fine=28, tol_rel=0.15
    )
    assert "dL_rel" in cert
