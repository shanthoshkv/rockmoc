"""Offline literature validation suite."""

from __future__ import annotations

from pathlib import Path

from rockmoc.validate import run_all_validation, validation_report_markdown


def test_all_validation_checks_pass():
    results = run_all_validation()
    failed = [r for r in results if not r.passed]
    assert not failed, "Failed checks:\n" + "\n".join(
        f"- {r.name}: {r.detail}" for r in failed
    )


def test_validation_report_writes(tmp_path: Path):
    md = validation_report_markdown()
    assert "NASA TM-X-1502" in md
    assert "PASS" in md or "FAIL" in md
    out = tmp_path / "validation_report.md"
    out.write_text(md, encoding="utf-8")
    assert out.exists()
