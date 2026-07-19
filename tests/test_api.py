"""Public API smoke tests."""

from __future__ import annotations

import rockmoc
from rockmoc import DesignSpec, PerfectGas, design_mln


def test_version():
    assert isinstance(rockmoc.__version__, str)
    assert rockmoc.__version__


def test_import_exports():
    assert callable(design_mln)
    assert PerfectGas is not None
    assert DesignSpec is not None


def test_quickstart_pattern():
    gas = PerfectGas(gamma=1.2, R=350.0)
    spec = DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=18)
    result = design_mln(spec)
    assert result.contour.length > 0
