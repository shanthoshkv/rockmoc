"""Smoke-test the deep validation runner (literature catalog)."""

from __future__ import annotations

from pathlib import Path

from rockmoc.validate.benchmark import run_deep_validation, write_validation_outputs

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "tests" / "data" / "validation_catalog.json"
OUT = ROOT / "tests" / "data" / "validation_results"


def test_deep_validation_catalog_exists():
    assert CATALOG.exists()


def test_deep_validation_all_regimes_pass():
    report = run_deep_validation(CATALOG)
    assert report["summary"]["n_cases"] >= 20
    failed = [c["case_id"] for c in report["cases"] if not c["passed"]]
    assert not failed, f"failed cases: {failed}"
    paths = write_validation_outputs(report, OUT)
    assert Path(paths["json"]).exists()
    assert Path(paths["csv"]).exists()
    assert Path(paths["markdown"]).exists()
