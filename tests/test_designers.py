"""End-to-end designer tests."""

from __future__ import annotations

import numpy as np
import pytest

from rockmoc import DesignSpec, PerfectGas, design_ideal, design_mln, design_tic
from rockmoc.exceptions import DesignError


@pytest.fixture
def gas():
    return PerfectGas(gamma=1.4, R=287.0)


def test_mln_basic(gas):
    spec = DesignSpec(
        gas=gas,
        throat_radius=0.05,
        exit_mach=2.4,
        n_characteristics=20,
        delta=1.0,
    )
    result = design_mln(spec)
    assert result.kind == "mln"
    assert result.contour.x.size >= 10
    assert result.metrics.length > 0
    assert result.metrics.exit_radius > spec.throat_radius
    assert result.metrics.area_ratio > 1.0
    assert np.all(np.isfinite(result.contour.x))
    assert np.all(np.isfinite(result.contour.r))
    assert np.all(result.contour.r > 0)
    # Geometric ε close to 1-D design ε with improved wall length model
    eps_1d = float(gas.area_ratio(2.4))
    assert abs(result.metrics.area_ratio - eps_1d) / eps_1d < 0.15
    L_over_Rt = result.metrics.length / spec.throat_radius
    assert 1.0 < L_over_Rt < 80.0
    assert np.all(np.diff(result.contour.r) >= -1e-9 * spec.throat_radius)


def test_mln_planar(gas):
    spec = DesignSpec(
        gas=gas,
        throat_radius=1.0,
        exit_mach=2.4,
        n_characteristics=25,
        delta=0.0,
    )
    result = design_mln(spec)
    assert "axisym_hf" in str(result.meta.get("method", "")) or result.meta.get(
        "designer"
    ) == "axisym_hf"
    assert result.metrics.length > 0


def test_ideal_basic(gas):
    spec = DesignSpec(
        gas=gas,
        throat_radius=0.04,
        exit_mach=3.0,
        n_characteristics=24,
        throat_rc_over_rt=1.5,
    )
    result = design_ideal(spec)
    assert result.kind == "ideal"
    assert result.metrics.length > 0
    mln = design_mln(
        DesignSpec(
            gas=gas,
            throat_radius=0.04,
            exit_mach=3.0,
            n_characteristics=24,
        )
    )
    assert result.metrics.length > 0 and mln.metrics.length > 0


def test_tic_shorter_than_ideal(gas):
    spec = DesignSpec(
        gas=gas,
        throat_radius=0.05,
        exit_mach=2.5,
        n_characteristics=20,
        truncate_fraction=0.7,
    )
    ideal = design_ideal(spec)
    tic = design_tic(spec)
    assert tic.kind == "tic"
    assert tic.metrics.length < ideal.metrics.length
    assert tic.meta.get("truncate_fraction") == 0.7


def test_area_ratio_targeting(gas):
    eps = float(gas.area_ratio(2.5))
    spec = DesignSpec(
        gas=gas,
        throat_radius=0.03,
        area_ratio=eps,
        n_characteristics=16,
    )
    result = design_mln(spec)
    assert abs(result.exit_mach_design - 2.5) < 1e-5


def test_invalid_spec():
    with pytest.raises(DesignError):
        DesignSpec(exit_mach=2.0, area_ratio=5.0).validate()
    with pytest.raises(DesignError):
        DesignSpec(throat_radius=-1, exit_mach=2.0).validate()


def test_export(tmp_path, gas):
    spec = DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.0, n_characteristics=12)
    result = design_mln(spec)
    csv_path = tmp_path / "wall.csv"
    json_path = tmp_path / "wall.json"
    result.export_csv(str(csv_path))
    result.export_json(str(json_path))
    assert csv_path.exists()
    assert json_path.exists()
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    assert data.shape[0] >= 2


def test_summary_keys(gas):
    result = design_mln(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.2, n_characteristics=15)
    )
    s = result.summary()
    assert "area_ratio" in s
    assert "Cf_vac_1d" in s
    assert s["kind"] == "mln"
