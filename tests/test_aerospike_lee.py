"""Tests for Lee / Vernacchia annular aerospike designer."""

from __future__ import annotations

import math

import numpy as np
import pytest

from rockmoc.aerospike import design_spike, design_spike_truncated, majdalani_mach_from_er
from rockmoc.thermo import PerfectGas


def test_majdalani_matches_area_mach():
    gas = PerfectGas(1.27)
    for er in [5.0, 12.0, 25.0]:
        Me_m = majdalani_mach_from_er(er, 1.27)
        Me_r = float(gas.mach_from_area_ratio(er, branch="supersonic"))
        assert abs(Me_m - Me_r) / Me_r < 0.02


def test_spike_default_vernacchia_params():
    """Default-ish params from mvernacc nozzle_solver.py."""
    r = design_spike(
        Pc=5e6,
        Tc=1900.0,
        gamma=1.27,
        molar_mass=20.18,
        Re=0.015,
        expansion_ratio=12.0,
        Pa=101325.0,
        n_points=2000,  # dense but faster than 10k in unit tests
    )
    s = r.summary()
    assert s["exit_mach"] > 1.0
    assert s["ht_m"] > 0.0
    assert s["At_m2"] > 0.0
    assert s["m_dot_kg_s"] > 0.0
    assert s["Isp_s"] > 0.0
    assert s["F_N"] > 0.0
    assert r.R[0] < r.Re  # spike starts inside lip
    assert r.X[-1] > r.X[0]
    assert np.all(np.isfinite(r.X))
    assert np.all(np.isfinite(r.R))
    # δ = π/2 − ν_e > 0
    assert 0.0 < r.delta_rad < 0.5 * math.pi
    # Radius decreases toward tip overall
    assert r.R[-1] < r.R[0]


def test_spike_truncation_reduces_length():
    # Vacuum ambient so cumulative Isp is monotonic increasing
    full = design_spike(
        Pc=5e6,
        Tc=1900.0,
        gamma=1.27,
        molar_mass=20.18,
        Re=0.015,
        expansion_ratio=12.0,
        Pa=0.0,
        n_points=1500,
    )
    trunc = design_spike_truncated(full, length_fraction=0.2)
    assert trunc.X[-1] < full.X[-1]
    assert trunc.Isp <= full.Isp + 1e-6
    assert trunc.meta["truncated"] is True


def test_dense_mesh_stable():
    # Keep ν_e < π/2 for the Lee/Vernacchia ideal model
    r = design_spike(
        Pc=3e6,
        Tc=3000.0,
        gamma=1.25,
        molar_mass=22.0,
        Re=0.05,
        expansion_ratio=20.0,
        Pa=1000.0,
        n_points=5000,
    )
    assert r.contour.x.size >= 1000
    assert r.Isp > 0


def test_high_expansion_rejected_when_nu_too_large():
    with pytest.raises(Exception):
        design_spike(
            Pc=3e6,
            Tc=3000.0,
            gamma=1.2,
            molar_mass=22.0,
            Re=0.05,
            expansion_ratio=80.0,
            Pa=1000.0,
            n_points=500,
        )
