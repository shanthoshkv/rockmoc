"""Free-shock separation: criterion formulas and physical limits."""

from __future__ import annotations

import pytest

from rockmoc.analysis.separation import psep_over_pa, separated_cf
from rockmoc.exceptions import DesignError
from rockmoc.thermo.perfect_gas import PerfectGas

GAS = PerfectGas(1.2)
CRITS = ("summerfield", "kalt_badal", "schilling", "schmucker", "stark", "morrisette")


def test_criterion_formulas():
    assert psep_over_pa("summerfield", M=3.0, p0_pa=50) == 0.4
    assert psep_over_pa("kalt_badal", M=3.0, p0_pa=32) == pytest.approx(2 / 3 / 2)
    assert psep_over_pa("schilling", M=3.0, p0_pa=1) == pytest.approx(0.582)
    assert psep_over_pa("schilling", M=3.0, p0_pa=1, conical=True) == pytest.approx(0.541)
    assert psep_over_pa("morrisette", M=3.0, p0_pa=50) == pytest.approx(0.3404, abs=1e-4)
    with pytest.raises(DesignError):
        psep_over_pa("nope", M=3.0, p0_pa=50)  # type: ignore[arg-type]


@pytest.mark.parametrize("crit", CRITS)
def test_limits(crit):
    kw = dict(area_ratio=40.0, pc=7e6, criterion=crit)
    vac = separated_cf(GAS, pa=0.0, **kw)
    assert not vac.separated and vac.Cf == vac.Cf_attached
    hi_alt = separated_cf(GAS, pa=1e3, **kw)  # p0/pa = 7000: attached
    assert not hi_alt.separated
    sl = separated_cf(GAS, pa=101325.0, **kw)  # ε=40 at sea level: separated
    assert sl.separated and 1.0 < sl.eps_sep < 40.0
    assert sl.Cf > sl.Cf_attached  # cutting off the over-expanded wall helps
    assert sl.psep_pa == pytest.approx(
        psep_over_pa(crit, M=sl.M_sep, p0_pa=7e6 / 101325.0), rel=1e-6)
    lower = separated_cf(GAS, pa=2e5, **kw)  # higher pa → separation moves upstream
    assert lower.eps_sep < sl.eps_sep


def test_separates_at_throat_raises():
    # p0/pa = 1: throat wall p/pa ≈ 0.56 < Kalt–Badal 2/3 → no supersonic flow
    with pytest.raises(DesignError):
        separated_cf(GAS, area_ratio=10.0, pc=1e5, pa=1e5, criterion="kalt_badal")
