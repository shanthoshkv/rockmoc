"""v0.7 features: kernels, lazy import, fallbacks, Rao MOC, convergence, exports, CLI, plots."""

from __future__ import annotations

import dataclasses
import math
import subprocess
import sys
import warnings

import numpy as np
import pytest

from rockmoc import DesignSpec, PerfectGas, design_ideal, design_mln
from rockmoc import _kernels as K
from rockmoc.analysis.convergence import richardson
from rockmoc.exceptions import MOCFallbackWarning
from rockmoc.moc.types import FlowPoint, PointKind


@pytest.fixture(autouse=True)
def _quiet_fallbacks():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MOCFallbackWarning)
        yield


def test_kernels_match_perfect_gas():
    for g in (1.15, 1.25, 1.4):
        gas = PerfectGas(g)
        for M in (1.0001, 1.3, 2.0, 4.5, 12.0):
            nu = float(gas.prandtl_meyer(np.array(M)))  # array path = reference
            assert K.pm(M, g) == pytest.approx(nu, rel=1e-12)
            assert K.pm_inv(nu, g) == pytest.approx(M, rel=1e-9)
            assert K.area_ratio(M, g) == pytest.approx(float(gas.area_ratio(np.array(M))), rel=1e-12)


def test_import_is_lazy():
    code = "import sys, rockmoc; print(any(m.startswith('rockmoc.designers') for m in sys.modules))"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"


def test_flowpoint_is_frozen():
    p = FlowPoint(0.0, 1.0, 0.1, 0.2, 1.5, 0.7, PointKind.INTERIOR)
    with pytest.raises(dataclasses.FrozenInstanceError):
        p.x = 1.0  # type: ignore[misc]
    assert p.copy(x=2.0).x == 2.0


def test_fallbacks_are_recorded_in_meta():
    r = design_ideal(DesignSpec(PerfectGas(1.4), throat_radius=0.05, exit_mach=2.4, n_characteristics=20))
    assert any("uniform near-sonic start" in w for w in r.meta["warnings"])


def test_richardson_recovers_order():
    n = np.array([10.0, 20.0, 40.0])
    st = richardson(list(1.0 + 3.0 / n**2))
    assert st["order"] == pytest.approx(2.0, abs=1e-9)
    assert st["extrapolated"] == pytest.approx(1.0, abs=1e-12)


def test_c_minus_lines_descend():
    r = design_mln(DesignSpec(PerfectGas(1.4), throat_radius=0.05, exit_mach=3.0, n_characteristics=16))
    lines = [ln for ln in r.net.c_minus if len(ln) > 2]
    assert lines and all(ln[-1].r < ln[0].r for ln in lines)


def test_rao_moc_hits_length_and_area():
    from rockmoc.moc.rao import design_rao_moc

    gas = PerfectGas(1.2)
    c, net, meta = design_rao_moc(gas, throat_radius=1.0, exit_mach=3.5, length_fraction=0.85, n_chars=48)
    assert c.x[-1] == pytest.approx(0.85 * meta["L15_m"], rel=1e-3)
    assert c.r[-1] == pytest.approx(math.sqrt(float(gas.area_ratio(3.5))), rel=1e-3)
    assert meta["theta_n_deg"] > meta["theta_e_deg"] > 0.0
    M = np.array([p.M for p in net.wall])
    assert np.all(np.diff(M) > -1e-6)


def test_field_cf_orders_designs():
    from rockmoc import field_integrated_cf
    from rockmoc.designers.rao_control import design_rao_control

    gas = PerfectGas(1.2)
    spec = DesignSpec(gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=24)
    cf = {}
    for name, res in (("mln", design_mln(spec)), ("rao", design_rao_control(spec, length_fraction=0.9))):
        cf[name] = field_integrated_cf(res.contour, gas, throat_radius=0.05, exit_mach_design=3.0,
                                       net=res.net).Cf_vac
    assert cf["rao"] < cf["mln"]  # divergence loss
    assert cf["rao"] > 0.97 * cf["mln"]


def test_exports(tmp_path):
    from rockmoc.export import export_blockmesh, export_cad_points, export_csv, export_dxf

    r = design_mln(DesignSpec(PerfectGas(1.4), throat_radius=0.05, exit_mach=2.5, n_characteristics=16))
    csv = np.loadtxt(export_csv(r.contour, tmp_path / "w.csv", units="mm"), delimiter=",", skiprows=1)
    assert csv[0, 1] == pytest.approx(50.0)
    xyz = np.loadtxt(export_cad_points(r.contour, tmp_path / "w.txt"))
    assert xyz.shape[1] == 3 and np.all(xyz[:, 2] == 0.0)
    dxf = export_dxf(r.contour, tmp_path / "w.dxf", mirror=True).read_text().split("\n")
    assert dxf[dxf.index("$ACADVER") + 2] == "AC1009"
    assert dxf.count("POLYLINE") == 2 and dxf.count("SEQEND") == 2 and dxf[-2] == "EOF"
    bm = export_blockmesh(r.contour, tmp_path / "system" / "blockMeshDict").read_text()
    for patch in ("inlet", "outlet", "wall", "front", "back", "axis"):
        assert f"    {patch} " in bm
    assert bm.count("spline") == 2


def test_cli_and_toml(tmp_path):
    from rockmoc.cli import main

    out = tmp_path / "out"
    spec = tmp_path / "n.toml"
    spec.write_text(
        '[gas]\ngamma = 1.25\n[design]\ntype = "tic"\nthroat_radius = 0.04\nexit_mach = 2.8\n'
        f'n_characteristics = 16\n[export]\nformats = ["csv", "dxf", "json"]\nout_dir = "{out.as_posix()}"\n'
    )
    assert main(["run", str(spec)]) == 0
    assert {p.name for p in out.iterdir()} == {"nozzle.csv", "nozzle.dxf", "nozzle.json"}
    assert main(["design", "--type", "mln", "--Me", "2.2", "-n", "12", "-o", str(out), "--name", "m"]) == 0
    assert (out / "m.csv").exists()


def test_plots_smoke(tmp_path):
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    from rockmoc.analysis import grid_convergence
    from rockmoc.plot import (
        plot_contour_overlay,
        plot_grid_convergence,
        plot_moc_field,
        plot_wall_properties,
    )

    gas = PerfectGas(1.4)
    before = dict(mpl.rcParams)
    r = design_mln(DesignSpec(gas, throat_radius=0.05, exit_mach=2.5, n_characteristics=12))
    plot_moc_field(r.contour, r.net, show=False, save_path=str(tmp_path / "f"), theme="light")
    plot_contour_overlay({"mln": r.contour}, show=False)
    plot_wall_properties(r.contour, gas, show=False)
    st = grid_convergence(gas, throat_radius=1.0, exit_mach=2.5, n_base=8, levels=3)
    plot_grid_convergence(st, show=False)
    assert (tmp_path / "f.svg").exists()  # no suffix → vector
    assert mpl.rcParams["axes.facecolor"] == before["axes.facecolor"]  # rc_context restored


def _source_state(x: float, r: float, gas: PerfectGas) -> tuple[float, float]:
    """Exact conical source flow from the origin: θ = atan(r/x), A/A* = (R/R*)², R* = 1."""
    R = math.hypot(x, r)
    M = float(gas.mach_from_area_ratio(R * R, branch="supersonic"))
    return math.atan2(r, x), M


@pytest.mark.parametrize("h", [0.04, 0.02, 0.01])
def test_axisym_interior_point_vs_exact_source_flow(h):
    """External analytic benchmark: one axisymmetric unit process inside an exact
    conical source flow. Local error is O(h³) (second-order method): observed
    7.1e-5, 8.9e-6, 1.1e-6, 1.4e-7 for h = 0.08 … 0.01."""
    from rockmoc.moc.mesh import interior_point

    gas = PerfectGas(1.4)
    xa, ra, xb, rb = 3.0, 1.0 + h, 3.0, 1.0 - h
    pts = []
    for x, r in ((xa, ra), (xb, rb)):
        th, M = _source_state(x, r, gas)
        pts.append(FlowPoint.from_mach_theta(x, r, M, th, gas))
    p = interior_point(pts[0], pts[1], gas, delta=1.0, max_iter=60, tol=1e-13)
    th_e, M_e = _source_state(p.x, p.r, gas)
    err = abs(p.theta - th_e) + abs(p.nu - float(gas.prandtl_meyer(M_e)))
    assert err < 0.2 * h**3
