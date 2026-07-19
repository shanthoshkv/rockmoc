#!/usr/bin/env python3
"""Run deep multi-regime validation and store reference vs solver error tables.

Outputs under ``tests/data/validation_results/``:

* ``validation_results.json`` — full structured report
* ``validation_results.csv``  — flat metric table (reference, solver, error %)
* ``validation_report.md``    — human-readable tables
* ``sources.md``              — literature source index

Usage::

    py -3.12 tests/run_deep_validation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rockmoc.validate.benchmark import main

if __name__ == "__main__":
    raise SystemExit(
        main(
            [
                "--catalog",
                str(ROOT / "tests" / "data" / "validation_catalog.json"),
                "--out",
                str(ROOT / "tests" / "data" / "validation_results"),
            ]
        )
    )
