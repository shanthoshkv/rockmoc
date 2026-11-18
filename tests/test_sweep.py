"""Sweep engine, new CLI commands, performance(), Contour.offset/from_csv, compare."""

from __future__ import annotations

import csv

import numpy as np
import pytest

from rockmoc.cli import main
from rockmoc.exceptions import DesignError
from rockmoc.geometry.contour import Contour
from rockmoc.sweep import parse_axis, run_sweep, write_rows

BASE = {
    "gas": {"gamma": 1.25, "R": 350.0},
    "design": {"type": "mln", "throat_radius": 0.05, "exit_mach": 2.5, "n_characteristics": 10},
    "performance": {"pc": 5e6, "Tc": 3200.0},
}


def test_parse_axis():
    assert parse_axis("Me=2:4:3") == ("Me", [2.0, 3.0, 4.0])
    assert parse_axis("type=mln,rao") == ("type", ["mln", "rao"])
    with pytest.raises(DesignError):
        parse_axis("Me")


def test_grid_pe_pc_error_row():
    rows = run_sweep(BASE, {"pc_bar": [40, 80], "pe": [2e4, 8e4]})
    assert len(rows) == 4 and not any("error" in r for r in rows)
    r = rows[0]  # pc=4e6, pe=2e4 → pe/pc = 0.005
    assert r["pc_Pa"] == pytest.approx(4e6)
    assert r["pe_pc"] == pytest.approx(0.005, rel=1e-6)
    # thrust = Cf * pc * At, mdot = pc*At/c*
    At = np.pi * 0.05**2
    assert r["thrust_vac_N"] == pytest.approx(r["Cf_vac"] * 4e6 * At)
    assert r["mdot_kg_s"] == pytest.approx(4e6 * At / r["cstar_m_s"])
    # pe >= pc → Me<=1 → recorded error, rest of grid survives
    bad = run_sweep(BASE, {"pe_pc": [0.01, 1.5]})
    assert "error" not in bad[0] and "error" in bad[1]
    assert write_rows(bad).splitlines()[0].startswith("case,pe_pc")


def test_sweep_rejects_two_exit_targets_and_unknown():
    with pytest.raises(DesignError):
        run_sweep(BASE, {"Me": [2, 3], "area_ratio": [10]})
    with pytest.raises(DesignError):
        run_sweep(BASE, {"bogus": [1]})


def test_of_sweep_uses_propellants():
    spec = {**BASE, "propellants": {"ox": "LOX", "fuel": "RP1", "of": 2.4}}
    rows = run_sweep(spec, {"of": [2.0, 3.0]})
    assert rows[0]["gamma_throat"] != rows[1]["gamma_throat"]
    assert rows[0]["gas_model"] in ("cea_equilibrium", "perfect")


def test_cli_sweep_validate_compare_stdout(tmp_path, capsys):
    out = tmp_path / "s.csv"
    main(["sweep", "-n", "10", "--vary", "Me=2:3:2", "--vary", "pc_bar=30,60", "--pc", "5e6",
          "--out", str(out), "--export", "csv", "--points", "20", "-o", str(tmp_path)])
    rows = list(csv.DictReader(out.open()))
    assert len(rows) == 4 and rows[0]["case"] == "1"
    wall = np.loadtxt(tmp_path / "case_0001.csv", delimiter=",", skiprows=1)
    assert len(wall) == 20

    toml = tmp_path / "n.toml"
    toml.write_text('[design]\nexit_mach = 2.5\nn_characteristics = 10\n')
    main(["validate", str(toml)])
    toml.write_text('[design]\nexit_mach = 2.5\narea_ratio = 9\n')
    with pytest.raises(DesignError):
        main(["validate", str(toml)])

    main(["compare", str(tmp_path / "case_0001.csv"), str(tmp_path / "case_0001.csv")])
    assert '"max_abs_dr_m": 0.0' in capsys.readouterr().out

    main(["design", "-n", "10", "--Me", "2.5", "--stdout", "-o", str(tmp_path)])
    assert capsys.readouterr().out.startswith("x_m,r_m")


def test_contour_offset_and_from_csv(tmp_path):
    x = np.linspace(0, 1, 40)
    c = Contour(x=x, r=0.1 + 0.0 * x)  # straight wall → offset is pure +r
    o = c.offset(0.01)
    assert np.allclose(o.r, 0.11) and np.allclose(o.x, x)
    p = tmp_path / "w.txt"
    rows = (f"{a * 1000} {b * 1000}" for a, b in zip(x, c.r, strict=True))
    p.write_text("# x_mm r_mm\n" + "\n".join(rows))
    back = Contour.from_csv(p, units="mm")
    assert np.allclose(back.x, x) and np.allclose(back.r, c.r)


def test_design_result_performance():
    from rockmoc.designers import DesignSpec, design_mln
    from rockmoc.thermo import PerfectGas

    r = design_mln(DesignSpec(PerfectGas(1.25, 350.0), exit_mach=2.5, n_characteristics=10))
    p = r.performance(pc=5e6, Tc=3200.0, pa=1.01e5)
    assert p["thrust_amb_N"] < p["thrust_vac_N"] and p["mdot_kg_s"] > 0
