"""Physics regression checks for MOC wall schedule, scaling, and performance."""

from __future__ import annotations

import math

import numpy as np
import pytest

from rockmoc.designers import DesignSpec, design_mln
from rockmoc.moc.axisym_hf import design_axisym_mln_hf
from rockmoc.performance.field_cf import field_integrated_cf
from rockmoc.performance.integrated import (
    divergence_efficiency,
    exit_kinetic_efficiency,
    integrated_performance,
)
from rockmoc.thermo.perfect_gas import PerfectGas


def test_canceling_wall_planar_invariant():
    """On planar canceling wall: θ + ν ≈ ν_e (classic MLN)."""
    gas = PerfectGas(1.4)
    Me = 3.0
    nu_e = float(gas.prandtl_meyer(Me))
    r = design_axisym_mln_hf(
        gas,
        throat_radius=0.05,
        exit_mach=Me,
        n_chars=24,
        delta=0.0,
        store_net=True,
    )
    assert r.net is not None and r.net.wall
    for p in r.net.wall:
        assert abs((p.theta + p.nu) - nu_e) < math.radians(1.5)


def test_canceling_wall_angle_matches_geometry():
    """Isotropic scale preserves wall angle: geometric θ ≈ stored θ."""
    gas = PerfectGas(1.25)
    r = design_axisym_mln_hf(
        gas, throat_radius=0.05, exit_mach=3.0, n_chars=32, store_net=False
    )
    c = r.contour
    th_g = np.arctan2(np.gradient(c.r), np.gradient(c.x))
    err = float(np.max(np.abs(th_g[2:-2] - np.asarray(c.theta)[2:-2])))
    assert err < math.radians(1.5)


def test_theta_max_is_half_nu_e_order():
    """Corner angle stays near ν_e/2 (not drifted by bogus axis-Me matching)."""
    gas = PerfectGas(1.4)
    nu_e = float(gas.prandtl_meyer(3.0))
    for n in (20, 48):
        r = design_axisym_mln_hf(
            gas, throat_radius=0.05, exit_mach=3.0, n_chars=n, store_net=False
        )
        # axisym mild bump to 0.52 ν_e
        assert 0.48 * nu_e < r.theta_max_rad < 0.55 * nu_e


def test_mln_exit_axial_theta_near_zero():
    gas = PerfectGas(1.25)
    res = design_mln(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=32),
        store_net=False,
    )
    assert abs(float(res.contour.theta[-1])) < math.radians(0.5)


def test_mln_length_mesh_trend():
    """Finer mesh should not collapse L; mild increase/stabilization expected."""
    gas = PerfectGas(1.4)
    Ls = []
    for n in (16, 32, 48):
        r = design_axisym_mln_hf(
            gas, throat_radius=0.05, exit_mach=3.0, n_chars=n, store_net=False
        )
        Ls.append(r.length / 0.05)
    assert Ls[0] > 2.0
    assert Ls[-1] > Ls[0] * 0.95  # not collapsing with n


def test_performance_no_false_wall_history_loss():
    """MLN with θ_e≈0 must have η_div=η_ke=1 (no mid-wall RMS penalty)."""
    gas = PerfectGas(1.4)
    res = design_mln(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=28),
        store_net=False,
    )
    assert divergence_efficiency(res.contour) > 0.999
    assert exit_kinetic_efficiency(res.contour) > 0.999
    perf = integrated_performance(
        res.contour, gas, throat_radius=0.05, exit_mach=3.0
    )
    assert perf.nozzle_efficiency > 0.999
    # Cf should match 1-D isentropic closely
    pe_pc = float(gas.pe_pc_from_mach(3.0))
    Cf_1d = gas.thrust_coefficient(
        pe_pc, 0.0, exit_mach=3.0, area_ratio=res.metrics.area_ratio
    )
    assert abs(perf.Cf_vac - Cf_1d) / Cf_1d < 0.01


def test_field_cf_reasonable():
    gas = PerfectGas(1.25, R=350.0)
    res = design_mln(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=28),
        store_net=True,
    )
    field = field_integrated_cf(
        res.contour,
        gas,
        throat_radius=0.05,
        exit_mach_design=3.0,
        pc=5e6,
        Tc=3200.0,
        net=res.net,
    )
    assert field.Cf_vac > 1.2
    assert field.Cf_vac < 2.2
    assert field.mass_flow > 0.0
    assert field.method == "contour_exit_station"


def test_area_ratio_matches_1d():
    gas = PerfectGas(1.4)
    r = design_axisym_mln_hf(
        gas, throat_radius=0.1, exit_mach=2.5, n_chars=24, store_net=False
    )
    eps_1d = float(gas.area_ratio(2.5))
    assert abs(r.area_ratio_geometric - eps_1d) / eps_1d < 1e-6
