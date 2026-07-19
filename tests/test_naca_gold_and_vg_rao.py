"""NACA gold regression, variable-γ along nozzle, Rao control-surface TOC."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from rockmoc import (
    DesignSpec,
    PerfectGas,
    design_mln,
    design_planar_mln_nasa,
)
from rockmoc.designers.rao_control import design_rao_control, rao_exit_angle
from rockmoc.thermo import VariableGammaGas, build_multi_gamma_table

DATA = Path(__file__).parent / "data" / "naca_rm_e8j12_gold.json"


def test_naca_rm_e8j12_gold_file():
    gold = json.loads(DATA.read_text(encoding="utf-8"))
    assert gold["gamma"] == 1.4
    gas = PerfectGas(1.4)
    for case in gold["cases"]:
        Me = case["Me"]
        nu = float(gas.prandtl_meyer(Me))
        assert abs(math.degrees(nu) - case["nu_e_deg"]) < 0.05
        assert abs(math.degrees(0.5 * nu) - case["theta_max_deg"]) < 0.05
        assert abs(float(gas.area_ratio(Me)) - case["area_ratio"]) < 0.01
        # Fine planar design
        r = design_planar_mln_nasa(Me, 1.4, dnu_deg=0.5, throat_half_height=1.0)
        assert abs(r.y_unit[-1] - case["area_ratio"]) < 0.02
        lo, hi = case["L_over_ht_range"]
        assert lo <= r.x_unit[-1] <= hi
        assert abs(r.theta_max_rad - 0.5 * nu) < 1e-9


def test_nasa_sample_in_gold():
    gold = json.loads(DATA.read_text(encoding="utf-8"))
    s = gold["nasa_tm_x_1502_sample"]
    r = design_planar_mln_nasa(s["Me"], s["gamma"], dnu_deg=0.5, throat_half_height=1.0)
    assert abs(r.y_unit[-1] - s["area_ratio"]) < 0.05
    assert r.x_unit[-1] > 0.5


def test_variable_gamma_along_nozzle_moc():
    """MOC with multi-γ table: γ_local varies with M along wall."""
    tab = build_multi_gamma_table(
        gamma_c=1.28, gamma_t=1.24, gamma_e=1.18, R=360.0, M_max=6.0, n=120
    )
    vg = VariableGammaGas(tab)
    # γ decreases along expansion
    assert vg.gamma_at_mach(1.1) > vg.gamma_at_mach(4.0)
    r = design_mln(
        DesignSpec(
            gas=vg,  # type: ignore[arg-type]
            throat_radius=0.05,
            exit_mach=3.0,
            n_characteristics=24,
        ),
        store_net=True,
    )
    assert r.metrics.length > 0
    # Wall γ(x) from table
    gx = vg.gamma_along_length(r.contour.x, x_throat=r.contour.x[0], x_exit=r.contour.x[-1])
    assert float(gx[0]) >= float(gx[-1]) - 1e-6  # non-increasing along nozzle
    # Net points should carry gamma_local when set
    if r.net and r.net.points:
        gl = [p.gamma_local for p in r.net.points if p.gamma_local == p.gamma_local]
        if gl:
            assert min(gl) < max(gl) or abs(min(gl) - max(gl)) < 0.15


def test_rao_control_surface():
    gas = PerfectGas(1.2)
    th = rao_exit_angle(3.5, 1.2, pa_pe=0.0)
    assert 0.0 < th < math.radians(25)
    r = design_rao_control(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.5, n_characteristics=40),
        length_fraction=0.8,
    )
    assert r.kind == "rao"
    assert r.meta["method"] == "rao_control_surface_toc"
    assert r.metrics.length > 0
    assert r.meta["theta_n_deg"] > r.meta["theta_e_rao_deg"]


def test_variable_gamma_rao_control():
    tab = build_multi_gamma_table(
        gamma_c=1.25, gamma_t=1.22, gamma_e=1.18, R=350.0, M_max=5.0, n=80
    )
    vg = VariableGammaGas(tab)
    r = design_rao_control(
        DesignSpec(gas=vg, throat_radius=0.04, exit_mach=3.0, n_characteristics=30),  # type: ignore[arg-type]
        length_fraction=0.75,
    )
    assert r.metrics.area_ratio > 1.0
