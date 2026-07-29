"""Tests for student-facing helpers: gases, conical, reports, tables, errors."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from rockmoc import (
    DesignSpec,
    compare_nozzle_designs,
    comparison_to_markdown,
    design_conical,
    design_mln,
    example_gas,
    explainer,
    generate_lab_report,
    list_example_gases,
)
from rockmoc.exceptions import DesignError
from rockmoc.student.limits import KNOWN_LIMITS


def test_example_gases():
    names = list_example_gases()
    assert "air" in names and "lox_rp1" in names
    g = example_gas("air")
    assert abs(g.gamma - 1.4) < 1e-9
    g2 = example_gas("lox_rp1")
    assert 1.15 < g2.gamma < 1.35


def test_design_spec_student():
    spec = DesignSpec.student(Me=3.0, gamma=1.25, Rt_mm=50, n=24)
    assert abs(spec.throat_radius - 0.05) < 1e-12
    assert spec.exit_mach == 3.0
    assert spec.n_characteristics == 24


def test_friendly_design_error_two_targets():
    with pytest.raises(DesignError, match="Fix:"):
        DesignSpec(
            gas=example_gas("air"),
            throat_radius=0.05,
            exit_mach=3.0,
            area_ratio=10.0,
        ).validate()


def test_conical_baseline():
    spec = DesignSpec.student(Me=3.0, gamma=1.4, Rt_mm=50, n=20)
    cone = design_conical(spec, half_angle_deg=15.0)
    assert cone.kind == "conical"
    assert cone.metrics.length > 0
    eps = cone.metrics.area_ratio
    assert abs(eps - float(spec.gas.area_ratio(3.0))) / eps < 1e-6
    # constant wall angle
    assert abs(math.degrees(cone.contour.theta[0]) - 15.0) < 1e-6


def test_report_and_lab_report(tmp_path: Path):
    spec = DesignSpec.student(Me=2.5, gamma=1.25, Rt_mm=40, n=20)
    mln = design_mln(spec, store_net=False)
    text = mln.report(length_unit="mm", angle_unit="deg")
    assert "rockmoc design report" in text
    assert "Limits:" in text
    assert "MLN" in text.upper() or "Minimum" in text
    assert "limitations" in mln.summary()

    path = mln.lab_report(
        tmp_path / "lab.md",
        student_name="Test Student",
        course="AE-101",
    )
    body = Path(path).read_text(encoding="utf-8")
    assert "Student checklist" in body
    assert "Test Student" in body
    assert "Known limitations" in body


def test_explainer():
    s = explainer("conical")
    assert "15" in s or "cone" in s.lower() or "Conical" in s


def test_comparison_exports(tmp_path: Path):
    spec = DesignSpec.student(Me=3.0, gamma=1.25, Rt_mm=50, n=20)
    cmp = compare_nozzle_designs(spec, include_conical=True)
    kinds = {r["kind"] for r in cmp.rows}
    assert "mln" in kinds and "conical" in kinds
    md = cmp.to_markdown()
    assert "L/Rt" in md or "ε" in md
    md2 = comparison_to_markdown(cmp)
    assert "|" in md2
    csv_path = tmp_path / "c.csv"
    cmp.to_csv(csv_path)
    assert csv_path.exists()
    tex = cmp.to_latex()
    assert "tabular" in tex


def test_known_limits_nonempty():
    assert len(KNOWN_LIMITS) >= 5
