"""v0.8 accuracy upgrades: integral BL, friction loss stack, Hall throat, multi-γ isentrope, Rao DE."""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

from rockmoc.analysis import boundary_layer as blmod
from rockmoc.analysis.boundary_layer import bartz_viscosity, correct_wall_bl
from rockmoc.designers import DesignSpec, design_mln
from rockmoc.exceptions import MOCFallbackWarning
from rockmoc.geometry.contour import Contour
from rockmoc.moc.rao import design_rao_moc
from rockmoc.moc.transonic import hall_uv, sauer_offset, sauer_uv
from rockmoc.performance.integrated import integrated_performance
from rockmoc.thermo.expansion_table import build_multi_gamma_table
from rockmoc.thermo.perfect_gas import PerfectGas
from rockmoc.thermo.variable_gamma import VariableGammaGas

GAS = PerfectGas(1.2, 350.0)


@pytest.fixture(scope="module")
def mln():
    return design_mln(DesignSpec(gas=GAS, throat_radius=0.05, exit_mach=3.5, n_characteristics=40),
                      store_net=False)


def _flat(M: float = 3.0) -> Contour:
    x = np.linspace(0.0, 1.0, 200)
    return Contour(x=x, r=np.full_like(x, 0.1), theta=np.zeros_like(x), M=np.full_like(x, M))


@pytest.mark.parametrize("turbulent", [True, False])
def test_integral_reduces_to_flat_plate_at_zero_pressure_gradient(turbulent):
    a = correct_wall_bl(_flat(), GAS, throat_radius=0.05, turbulent=turbulent, method="integral")
    b = correct_wall_bl(_flat(), GAS, throat_radius=0.05, turbulent=turbulent, method="empirical")
    assert a.meta["method"] == "momentum_integral"
    assert b.meta["method"] == "empirical_flat_plate"
    np.testing.assert_allclose(a.theta_mom, b.theta_mom, rtol=1e-4)


def test_integral_failure_falls_back_to_correlation(mln, monkeypatch):
    def boom(*a, **k):
        raise ValueError("forced")

    monkeypatch.setattr(blmod, "solve_ivp", boom)
    with pytest.warns(MOCFallbackWarning, match="flat-plate correlation used"):
        bl = correct_wall_bl(mln.contour, GAS, throat_radius=0.05)
    assert bl.meta["method"] == "empirical_flat_plate"
    assert bl.meta["fallback_reason"] == "forced"
    assert bl.delta_star[-1] > 0


def test_cold_wall_thins_displacement_thickness(mln):
    hot = correct_wall_bl(mln.contour, GAS, throat_radius=0.05, pc=7e6, Tc=3400.0)
    cold = correct_wall_bl(mln.contour, GAS, throat_radius=0.05, pc=7e6, Tc=3400.0, wall_temperature=800.0)
    assert cold.delta_star[-1] < hot.delta_star[-1]
    assert cold.H[-1] < hot.H[-1]


def test_friction_loss_is_realistic_and_enters_cf(mln):
    bl = correct_wall_bl(mln.contour, GAS, throat_radius=0.05, pc=7e6, Tc=3400.0)
    kw = dict(throat_radius=0.05, exit_mach=3.5, pc=7e6, Tc=3400.0)
    clean = integrated_performance(mln.contour, GAS, **kw)
    lossy = integrated_performance(mln.contour, GAS, **kw, delta_cf_friction=bl.delta_cf_friction,
                                   eta_chemical=0.98)
    frac = bl.delta_cf_friction / clean.Cf_vac
    assert 0.002 < frac < 0.03  # 0.2–3 % for a 5 cm throat
    assert lossy.Cf_vac == pytest.approx(clean.Cf_vac * 0.98 - bl.delta_cf_friction, rel=1e-12)
    assert lossy.meta["loss_stack"]["eta_chemical"] == 0.98


@pytest.mark.parametrize("M", [2.0, 3.0, 5.0])
@pytest.mark.parametrize("tw_over_taw", [None, 0.5])
def test_skin_friction_matches_van_driest_ii(M, tw_over_taw):
    """Eckert reference-T closure vs Van Driest II (same incompressible law, μ ∝ T^0.6)."""
    g, Tc, pc = GAS.gamma, 3400.0, 7e6
    m = 0.5 * (g - 1.0) * M * M
    Te = Tc / (1.0 + m)
    rec = (4.0 * g / (9.0 * g - 5.0)) ** (1.0 / 3.0)
    Taw = Te * (1.0 + rec * m)
    Tw = None if tw_over_taw is None else tw_over_taw * Taw
    bl = correct_wall_bl(_flat(M), GAS, throat_radius=0.05, pc=pc, Tc=Tc, wall_temperature=Tw)
    i = len(bl.theta_mom) // 2
    pe = pc / (1.0 + m) ** (g / (g - 1.0))
    rho_e = pe / (GAS.R * Te)
    U = M * math.sqrt(g * GAS.R * Te)
    Re_th = rho_e * U * bl.theta_mom[i] / float(bartz_viscosity(Te, 8314.462618 / GAS.R))
    tw = (Taw if Tw is None else Tw) / Te
    taw = Taw / Te
    a = math.sqrt(rec * m / tw)
    b = taw / tw - 1.0
    d = math.sqrt(b * b + 4.0 * a * a)
    Fc = (taw - 1.0) / (math.asin((2.0 * a * a - b) / d) + math.asin(b / d)) ** 2
    cf_vd = 0.0256 * (Re_th / tw**0.6) ** -0.25 / Fc
    assert bl.cf[i] == pytest.approx(cf_vd, rel=0.03)


def test_bartz_viscosity_magnitude():
    assert bartz_viscosity(3500.0, 22.0) == pytest.approx(7.45e-5, rel=0.01)


def _potential_residual(uv, R: float, g: float = 1.2) -> float:
    """Max |continuity residual| of the exact axisymmetric potential equation."""
    h, worst = 1e-6, 0.0
    for x in np.linspace(-0.3, 0.3, 5) / math.sqrt(R):
        for r in np.linspace(0.1, 0.9, 5):
            u, v = uv(x, r, R)
            ux = (uv(x + h, r, R)[0] - uv(x - h, r, R)[0]) / (2 * h)
            vx = (uv(x + h, r, R)[1] - uv(x - h, r, R)[1]) / (2 * h)
            ur = (uv(x, r + h, R)[0] - uv(x, r - h, R)[0]) / (2 * h)
            vr = (uv(x, r + h, R)[1] - uv(x, r - h, R)[1]) / (2 * h)
            a2 = 0.5 * (g + 1) - 0.5 * (g - 1) * (u * u + v * v)
            worst = max(worst, abs((a2 - u * u) * ux - u * v * (ur + vx) + (a2 - v * v) * vr + a2 * v / r))
    return worst


def test_hall_converges_faster_than_sauer():
    g = 1.2
    hall = lambda x, r, R: hall_uv(x, r, g, 1.0, R)  # noqa: E731

    def sauer(x, r, R):
        u, v = sauer_uv(x - sauer_offset(g, 1.0, R), r, g, 1.0, R)
        return 1.0 + u, v

    for R in (4.0, 8.0):
        assert _potential_residual(hall, R) < 0.5 * _potential_residual(sauer, R)
    # observed order: Hall ~ε^3.4, Sauer ~ε^2.5
    p_hall = math.log2(_potential_residual(hall, 8.0) / _potential_residual(hall, 16.0))
    assert p_hall > 3.0


def test_multi_gamma_constant_reduces_to_perfect_gas():
    t = build_multi_gamma_table(gamma_c=1.2, gamma_t=1.2, gamma_e=1.2, R=350.0, Tc=3000.0)
    np.testing.assert_allclose(t.T_Tc * GAS.temperature_ratio(t.M), 1.0, atol=1e-6)
    np.testing.assert_allclose(t.area_ratio / GAS.area_ratio(t.M), 1.0, atol=1e-6)
    np.testing.assert_allclose(t.nu, GAS.prandtl_meyer(t.M), atol=1e-5)
    vg = VariableGammaGas(t)
    assert vg.characteristic_velocity(3300.0) == pytest.approx(GAS.characteristic_velocity(3300.0), rel=1e-6)


def test_rao_coarse_lattice_reaches_target_length():
    # n = 20 used to raise "not reachable": DE crossed θ = μ (acos branch singularity)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MOCFallbackWarning)
        _, _, coarse = design_rao_moc(GAS, throat_radius=0.05, exit_mach=3.5, n_chars=20)
        _, _, fine = design_rao_moc(GAS, throat_radius=0.05, exit_mach=3.5, n_chars=80)
    assert coarse["length_m"] / coarse["L15_m"] == pytest.approx(0.80, abs=1e-3)
    assert abs(coarse["theta_n_deg"] - fine["theta_n_deg"]) < 0.6
    assert abs(coarse["theta_e_deg"] - fine["theta_e_deg"]) < 0.5
