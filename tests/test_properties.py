"""Property-based tests (Hypothesis): unit-process invariants and MLN design properties."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from rockmoc import PerfectGas
from rockmoc.exceptions import MOCFallbackWarning
from rockmoc.moc.types import FlowPoint

pytest.importorskip("hypothesis", reason="install rockmoc[dev] for property tests")
from hypothesis import given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402


@pytest.fixture(autouse=True)
def _quiet_fallbacks():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MOCFallbackWarning)
        yield


@settings(max_examples=60, deadline=None)
@given(
    g=st.floats(1.15, 1.4),
    Ma=st.floats(1.3, 3.0),
    Mb=st.floats(1.3, 3.0),
    tha=st.floats(-0.2, 0.3),
    thb=st.floats(-0.2, 0.3),
)
def test_planar_unit_process_invariants(g, Ma, Mb, tha, thb):
    """δ = 0: (θ+ν) is carried along C− from A and (θ−ν) along C+ from B, exactly."""
    from rockmoc.exceptions import MeshError, ThermoError
    from rockmoc.moc.mesh import interior_point

    gas = PerfectGas(g)
    a = FlowPoint.from_mach_theta(0.0, 1.0, Ma, tha, gas)
    b = FlowPoint.from_mach_theta(0.05, 0.5, Mb, thb, gas)
    try:
        p = interior_point(a, b, gas, delta=0.0, max_iter=60, tol=1e-12)
    except (MeshError, ThermoError):
        return
    if (a.theta - a.mu) - (b.theta + b.mu) >= -1e-3:  # characteristics do not intersect downstream
        return
    assert p.theta + p.nu == pytest.approx(a.theta + a.nu, abs=1e-9)
    assert p.theta - p.nu == pytest.approx(b.theta - b.nu, abs=1e-9)


@settings(max_examples=6, deadline=None)
@given(g=st.floats(1.18, 1.4), Me=st.floats(1.8, 3.5))
def test_axisym_mln_properties(g, Me):
    """Axis reaches Me, exit area is mass-exact, exit is axial, wall Mach rises."""
    from rockmoc.moc.axisym_hf import design_axisym_mln_hf

    r = design_axisym_mln_hf(PerfectGas(g), throat_radius=1.0, exit_mach=Me, n_chars=16,
                             match_exit_area=False)
    assert r.axis_exit_mach == pytest.approx(Me, rel=1e-6)
    assert r.area_ratio_geometric == pytest.approx(r.area_ratio_1d, rel=1e-6)
    assert abs(float(r.contour.theta[-1])) < 1e-12
    M = np.array([p.M for p in r.net.wall[1:]])
    assert np.all(np.diff(M) > -1e-6)
    assert r.theta_max_rad < 0.5 * float(PerfectGas(g).prandtl_meyer(Me))
