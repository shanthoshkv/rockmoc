"""Tests for features 1–8 (BL, tables, Rao, throat, performance, spike extras)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from rockmoc import (
    DesignSpec,
    PerfectGas,
    compare_nozzle_designs,
    correct_wall_bl,
    design_mln,
    design_rao,
    design_spike,
    design_spike_with_cowl,
    integrated_performance,
    off_design_altitude_sweep,
    sauer_throat_ivl,
)
from rockmoc.aerospike import truncated_spike_with_base
from rockmoc.thermo.expansion_table import build_multi_gamma_table, build_perfect_gas_table


def test_bl_correction():
    gas = PerfectGas(1.4)
    r = design_mln(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.5, n_characteristics=20)
    )
    bl = correct_wall_bl(r.contour, gas, throat_radius=0.05, pc=4e6, Tc=2800.0)
    assert bl.delta_star[-1] > 0
    assert bl.corrected.r[-1] >= bl.inviscid.r[-1] - 1e-12
    assert "delta_star_exit_m" in bl.summary()


def test_expansion_tables():
    gas = PerfectGas(1.3)
    t = build_perfect_gas_table(gas, M_max=4.0, n=50)
    assert t.state_at_mach(2.0)["area_ratio"] > 1.0
    m = build_multi_gamma_table(
        gamma_c=1.25, gamma_t=1.22, gamma_e=1.18, R=360.0, M_max=4.0, n=40
    )
    assert m.mode == "multi_gamma"
    assert m.mach_from_nu(float(m.nu[10])) > 1.0


def test_sauer_throat():
    gas = PerfectGas(1.25)
    ivl = sauer_throat_ivl(gas, throat_radius=0.04, rc_over_rt=2.0, n_points=15)
    assert np.all(ivl.M > 1.0)
    assert ivl.points[0].r == pytest.approx(0.0, abs=1e-12)


def test_rao():
    gas = PerfectGas(1.2)
    r = design_rao(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.5, n_characteristics=40),
        length_fraction=0.8,
    )
    assert r.kind == "rao"
    assert r.metrics.length > 0
    assert r.meta["theta_n_deg"] > r.meta["theta_e_deg"]
    assert abs(r.meta["length_fraction"] - 0.8) < 1e-9


def test_integrated_and_offdesign():
    gas = PerfectGas(1.25)
    r = design_mln(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=20)
    )
    perf = integrated_performance(
        r.contour, gas, throat_radius=0.05, exit_mach=3.0, pc=6e6, Tc=3100.0, pa=101325.0
    )
    assert perf.Cf_vac > 1.0
    assert 0.9 < perf.nozzle_efficiency <= 1.0
    rows = off_design_altitude_sweep(
        r.contour, gas, throat_radius=0.05, exit_mach=3.0, altitudes_m=[0, 10e3, 40e3]
    )
    assert len(rows) == 3


def test_compare_designs():
    gas = PerfectGas(1.25)
    c = compare_nozzle_designs(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.5, n_characteristics=18)
    )
    assert len(c.rows) >= 4  # mln, ideal, tic, rao (+ conical by default)
    kinds = {r["kind"] for r in c.rows}
    assert {"mln", "ideal", "tic", "rao"}.issubset(kinds)
    assert "rao" in c.summary_table()


def test_spike_base_and_cowl():
    full = design_spike(
        Pc=5e6, Tc=1900, gamma=1.27, molar_mass=20.18,
        Re=0.015, expansion_ratio=12, Pa=5e4, n_points=500,
    )
    trunc, bp, perf = truncated_spike_with_base(full, length_fraction=0.2)
    assert bp.p_base >= 0
    assert math.isfinite(perf["F_corrected_N"])
    stack = design_spike_with_cowl(
        Pc=5e6, Tc=1900, gamma=1.27, molar_mass=20.18,
        Re=0.015, expansion_ratio=12, n_points=500, primary_eps=1.8, truncate=0.25,
    )
    assert stack["cowl"].area_ratio_primary == 1.8
    assert "base_performance" in stack
