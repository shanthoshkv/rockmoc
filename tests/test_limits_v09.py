"""v0.9: fixes for the documented limits (throat, MLN data, separation, BL, shocks, kinetics, two-phase)."""

from __future__ import annotations

import math
import warnings

import pytest

from rockmoc.analysis.boundary_layer import two_layer_thicknesses
from rockmoc.analysis.separation import psep_over_pa, separated_cf
from rockmoc.exceptions import MOCFallbackWarning
from rockmoc.moc.transonic import discharge_coefficient, kl_uv, mach_from_mstar
from rockmoc.thermo.perfect_gas import PerfectGas


def _M(uv, g=1.4):
    return mach_from_mstar(math.hypot(*uv), g)


def test_kliegel_levine_matches_cuffel_tight_throat():
    assert _M(kl_uv(0.0, 0.0, 1.4, 1.0, 0.625)) == pytest.approx(0.8, abs=0.03)
    assert discharge_coefficient(1.4, 0.625) == pytest.approx(0.985, abs=0.004)


def test_kliegel_levine_axial_scaling_is_irrotational_to_truncation():
    h, err = 1e-6, 0.0
    for x in (-0.1, 0.0, 0.1):
        for r in (0.2, 0.5, 0.8):
            ur = (kl_uv(x, r + h, 1.4, 1, 16)[0] - kl_uv(x, r - h, 1.4, 1, 16)[0]) / (2 * h)
            vx = (kl_uv(x + h, r, 1.4, 1, 16)[1] - kl_uv(x - h, r, 1.4, 1, 16)[1]) / (2 * h)
            err = max(err, abs(ur - vx))
    assert err < 1e-4  # O(ε⁴); the √(R+1) scaling would leave O(ε²) ≈ 1e-3


def test_two_layer_shape_factor_matches_bartz_magnitudes():
    assert two_layer_thicknesses(3.0, 0.5, 1.0, 1.4)[0] == pytest.approx(3.57, abs=0.02)
    assert two_layer_thicknesses(3.0, 0.5, 2.0, 1.4)[0] == pytest.approx(1.40, abs=0.02)
    assert two_layer_thicknesses(3.0, 0.5, 4.0, 1.4)[0] < 0.0  # Bartz: δ*/θ can go negative


def test_separation_criteria_and_rocketcea_cross_check():
    assert psep_over_pa("stark", M=3.0, p0_pa=50) == pytest.approx(math.pi / 9)
    assert psep_over_pa("schmucker", M=3.0, p0_pa=50) == pytest.approx(4.64 ** -0.64)
    rc = pytest.importorskip("rocketcea.separated_Cf")
    ref = rc.ambientCf(gam=1.2, epsTot=40.0, Pc=1000.0, Pamb=14.7)[0]
    r = separated_cf(PerfectGas(1.2), area_ratio=40.0, pc=1000.0, pa=14.7, criterion="kalt_badal", recovery="kalt_badal")
    assert r.separated and r.Cf == pytest.approx(ref, rel=1e-6)
    att = separated_cf(PerfectGas(1.2), area_ratio=10.0, pc=7e6, pa=0.0)
    assert not att.separated


def test_internal_shock_detector_flags_top_not_ideal():
    from rockmoc.designers import DesignSpec, design_ideal, design_rao
    from rockmoc.moc.analysis import analyze_wall

    gas = PerfectGas(1.2, 350.0)
    spec = DesignSpec(gas=gas, throat_radius=1.0, exit_mach=3.5, n_characteristics=40)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MOCFallbackWarning)
        assert analyze_wall(design_rao(spec).contour, gas, rc_over_rt=1.5).has_internal_shock
        assert not analyze_wall(design_ideal(spec).contour, gas, rc_over_rt=1.5).has_internal_shock


def _bell(scale):
    from rockmoc.designers import DesignSpec, design_rao_control

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MOCFallbackWarning)
        c = design_rao_control(DesignSpec(gas=PerfectGas(1.2, 350.0), throat_radius=1.0, exit_mach=3.0,
                                          n_characteristics=40)).contour
    return (c.x - c.x[0]) * scale, c.r * scale


def test_two_phase_no_lag_limit_and_trend():
    from rockmoc.performance.two_phase import two_phase_nozzle

    x, r = _bell(0.05)
    kw = dict(x=x, r=r, gamma_gas=1.2, R_gas=350.0, pc=7e6, Tc=3400.0, particle_mass_fraction=0.3)
    assert two_phase_nozzle(diameter=1e-8, **kw).eta_two_phase == pytest.approx(1.0, abs=1e-4)
    e3, e10 = (two_phase_nozzle(diameter=d, **kw).eta_two_phase for d in (3e-6, 10e-6))
    assert 0.9 < e10 < e3 < 1.0


def test_kinetics_limits():
    pytest.importorskip("cantera")
    from rockmoc.thermo.kinetics import kinetic_nozzle

    x, r = _bell(0.05)
    k0 = kinetic_nozzle(reactants="H2:3.0, O2:1.0", pc=3e6, x=x, r=r, rate_multiplier=0.0)
    kf = kinetic_nozzle(reactants="H2:3.0, O2:1.0", pc=3e6, x=x, r=r, rate_multiplier=1e3)
    k = kinetic_nozzle(reactants="H2:3.0, O2:1.0", pc=3e6, x=x, r=r)
    assert k0.Isp_kinetic == pytest.approx(k0.Isp_frozen, rel=1e-4)
    assert kf.Isp_kinetic == pytest.approx(kf.Isp_shifting, rel=2e-3)
    assert k0.Isp_frozen < k.Isp_kinetic < k.Isp_shifting


@pytest.mark.slow
def test_euler_throat_matches_kliegel_levine_at_rc2():
    from rockmoc.moc.transonic_euler import cuffel_wall, solve_throat_extrapolated

    r = solve_throat_extrapolated(cuffel_wall(2.0, r_inlet=2.0), -2.2, 1.0)
    assert r["Cd"] == pytest.approx(discharge_coefficient(1.4, 2.0), abs=0.003)
    assert r["M_axis"] == pytest.approx(_M(kl_uv(0.0, 0.0, 1.4, 1.0, 2.0)), abs=0.02)
