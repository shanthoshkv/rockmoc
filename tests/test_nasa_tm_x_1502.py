"""Golden / regression tests for NASA TM-X-1502 planar MLN."""

from __future__ import annotations

import math

import numpy as np
import pytest

from rockmoc.moc.planar_nasa import design_planar_mln_nasa, mesh_convergence_planar
from rockmoc.thermo import PerfectGas


def test_sample_case_runs_me405_g136():
    """Sample input from TM-X-1502: ME=4.05, GAM=1.36 (Δν=0.5° for CI speed).

    Paper sample uses DV=0.10°; production should use ≤0.1°. Tests use 0.5°
    so the O(n²) characteristic net stays tractable.
    """
    r = design_planar_mln_nasa(4.05, 1.36, dnu_deg=0.5, throat_half_height=1.0)
    s = r.summary()
    assert s["method"].startswith("NASA_TM-X-1502")
    assert r.contour.x.size >= 10
    assert r.x_unit[0] == pytest.approx(0.0, abs=1e-12)
    assert r.y_unit[0] == pytest.approx(1.0, abs=1e-12)
    # Diverging: exit height > throat
    assert r.y_unit[-1] > 1.0
    assert r.x_unit[-1] > 0.0
    # θ_max = ν_e/2
    gas = PerfectGas(1.36)
    nu_e = float(gas.prandtl_meyer(4.05))
    assert abs(r.theta_max_rad - 0.5 * nu_e) < 1e-12
    # Final wall angle near axial
    assert abs(r.contour.theta[-1]) < math.radians(1.0)
    # Physical axial length (was near-zero before wall geometry fix)
    assert r.x_unit[-1] > 0.5
    assert r.summary()["length_m"] > 0.0


def test_exit_height_matches_area_ratio():
    """Uniform-exit planar design: y_e/y_t = A/A*(M_e) (TM-X-1502)."""
    gas = PerfectGas(1.4)
    Me = 2.4
    eps = float(gas.area_ratio(Me))
    r_f = design_planar_mln_nasa(Me, 1.4, dnu_deg=0.5)
    assert r_f.y_unit[-1] == pytest.approx(eps, rel=1e-9)
    assert np.isfinite(r_f.x_unit[-1])
    assert r_f.x_unit[-1] > 0.5  # physical MLN length / ht


def test_physical_scale():
    r = design_planar_mln_nasa(2.0, 1.4, dnu_deg=0.5, throat_half_height=0.05)
    assert r.contour.r[0] == pytest.approx(0.05, rel=1e-9)
    assert r.contour.x[0] == pytest.approx(0.0, abs=1e-12)


def test_mesh_convergence_monotonic_trend():
    rows = mesh_convergence_planar(3.0, 1.4, dnu_deg_list=[1.0, 0.5, 0.25])
    assert len(rows) == 3
    for row in rows:
        assert row["exit_x_over_ht"] > 0
        assert row["exit_y_over_ht"] > 1.0


def test_theta_max_rule():
    for Me, g in [(2.0, 1.4), (3.5, 1.25), (4.05, 1.36)]:
        r = design_planar_mln_nasa(Me, g, dnu_deg=0.5)
        gas = PerfectGas(g)
        assert abs(r.theta_max_rad - 0.5 * float(gas.prandtl_meyer(Me))) < 1e-10
