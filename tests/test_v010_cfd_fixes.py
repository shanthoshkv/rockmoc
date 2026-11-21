"""v0.10: fixes found by the OpenFOAM validation campaign (cfd_validation/FINDINGS.md)."""

import math

import numpy as np
import pytest

from rockmoc import DesignSpec, design_conical, design_rao, design_rao_toc
from rockmoc.analysis.separation import separated_cf
from rockmoc.thermo import PerfectGas

GAS = PerfectGas(1.22, 355.0)


def _spec():
    return DesignSpec(gas=GAS, throat_radius=0.05, exit_mach=3.0, n_characteristics=40)


def test_cone_divergence_counted_once():
    r = design_conical(_spec(), half_angle_deg=15.0)
    p = r.performance(pc=5e6, Tc=3200.0)
    lam = 0.5 * (1 + math.cos(math.radians(15.0)))
    assert p["geometric_model"] == "conical_source_flow"
    assert p["Cf_vac"] == pytest.approx(p["Cf_1d_vac"] * lam, rel=1e-9)  # CFD: within 0.3 %


def test_rao_exit_angle_grows_as_bell_shortens():
    th = [design_rao(_spec(), length_fraction=f).meta["theta_e_deg"] for f in (1.0, 0.8, 0.6)]
    assert th[0] < th[1] < th[2]
    c = design_rao(_spec(), length_fraction=0.8).contour
    slope = math.degrees(math.atan((c.r[-1] - c.r[-2]) / (c.x[-1] - c.x[-2])))
    assert math.degrees(c.theta[-1]) == pytest.approx(slope, abs=0.3)  # θ consistent with geometry


def test_short_bell_loses_more_than_long():
    eta = [
        design_rao(_spec(), length_fraction=f).performance()["loss_stack"]["eta_geometric"]
        for f in (1.0, 0.6)
    ]
    assert eta[1] < eta[0] - 0.02  # MOC wall integral sees the short-bell over-expansion


def test_toc_theta_matches_geometry_and_keeps_throat_arc():
    c = design_rao_toc(_spec()).contour
    dth = np.degrees(np.arctan(np.gradient(c.r, c.x)))[5:-5] - np.degrees(np.asarray(c.theta))[5:-5]
    assert np.max(np.abs(dth)) < 0.5
    k = (c.x > 0) & (c.x < 0.2 * 0.05)
    assert np.median(c.x[k] ** 2 / (2 * (c.r[k] - 0.05))) / 0.05 == pytest.approx(1.5, rel=0.05)


def test_separation_uses_given_wall_pressure():
    kw = dict(area_ratio=6.34, pc=1e6, pa=101325.0, criterion="schmucker")
    base = separated_cf(GAS, **kw)
    e = np.linspace(1.01, 6.34, 50)
    p1d = np.array(
        [float(GAS.pe_pc_from_mach(GAS.mach_from_area_ratio(x, branch="supersonic"))) for x in e]
    )
    M = np.array([float(GAS.mach_from_area_ratio(x, branch="supersonic")) for x in e])
    hi = separated_cf(
        GAS, wall_eps=e, wall_p_pc=1.2 * p1d, wall_M=M, **kw
    )  # wall above 1-D → separates later
    assert base.separated and hi.eps_sep > base.eps_sep
    assert base.Cf_vac_attached > base.Cf  # vacuum part excludes the ambient pressure term


def test_ideal_throat_arc_is_densely_sampled():
    from rockmoc import design_ideal

    c = design_ideal(_spec()).contour
    sl = np.degrees(np.arctan2(np.diff(c.r), np.diff(c.x)))
    assert np.max(np.abs(np.diff(sl))) < 2.0  # v0.9: 5.9° polyline kinks on the throat arc


def test_bl_corrected_loss_uses_exit_control_volume():
    from rockmoc import design_production

    s = DesignSpec(gas=GAS, throat_radius=0.05, exit_mach=3.0, n_characteristics=60)
    p = design_production(s, apply_bl=True, kind="mln", pc=1e6, Tc=3200.0, certificate=False)
    m = p.bl.meta
    assert 0.0 < m["delta_cf_exit_vac"] < m["delta_cf_friction"]  # p_e on displaced wall offsets part of τ_w
    assert p.field.summary()["Isp_vac_s"] == pytest.approx(276.45, rel=3e-3)  # OpenFOAM v2 y+~0.5


def test_uncorrected_wall_bl_loss_matches_cfd():
    from rockmoc import design_mln
    from rockmoc.analysis.boundary_layer import correct_wall_bl

    d = design_mln(DesignSpec(gas=GAS, throat_radius=0.05, exit_mach=3.0, n_characteristics=60))
    m = correct_wall_bl(d.contour, GAS, throat_radius=0.05, pc=1e6, Tc=3200.0, turbulent=True).meta
    isp = d.performance(pc=1e6, Tc=3200.0, delta_cf_friction=m["delta_cf_uncorrected_vac"])["Isp_vac_s"]
    assert isp == pytest.approx(275.55, rel=3e-3)  # OpenFOAM v2, y+~0.5, wall not BL-corrected


def test_planar_mln_wall_from_characteristics():
    from rockmoc import design_planar_mln_nasa

    p = design_planar_mln_nasa(2.4, 1.4)
    assert p.meta["area_ratio_moc_raw"] == pytest.approx(p.meta["area_ratio_1d"], rel=1e-3)  # mass, unforced
    assert p.contour.x[-1] == pytest.approx(8.087, rel=2e-3)  # independent Anderson MOC; v0.9 gave 8.715


def test_transonic_series_choice_matches_cfd_throat_plane():
    from rockmoc.moc.transonic import HALL_MIN_RC_RT, kl_uv

    g = 1.22
    M = lambda V: math.sqrt(2 * V * V / ((g + 1) - (g - 1) * V * V))  # noqa: E731
    assert HALL_MIN_RC_RT > 2.0  # KL used for the common Rc/Rt = 1-2 throats
    # OpenFOAM v2 throat-plane axis Mach: Rc/Rt 1.5 -> 0.880, 1 -> 0.865, 0.625 -> 0.814
    for R, cfd in ((1.5, 0.880), (1.0, 0.865), (0.625, 0.814)):
        assert M(math.hypot(*kl_uv(0.0, 0.0, g, 1.0, R))) == pytest.approx(cfd, abs=0.02)


def test_rao80_sea_level_cf_matches_cfd():
    r = design_rao(DesignSpec(gas=GAS, throat_radius=0.05, exit_mach=3.0, n_characteristics=40), length_fraction=0.8)
    p = r.performance(pc=1e6, Tc=3200.0, pa=101325.0, separation="schmucker")
    assert p["Cf_amb"] == pytest.approx(0.997, rel=0.01)  # OpenFOAM v2, y+~1: lip-only separation, Cf 0.997
