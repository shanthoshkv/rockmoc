"""Gated comparisons against published experiments / reference designs (docs/validation_experimental.md)."""

import pytest

from rockmoc.moc.transonic import hall_line
from rockmoc.validate.experimental import run_experimental_validation


@pytest.mark.slow
def test_published_data_gated_checks_pass():
    failed = [f"{c.case}: {c.quantity} ref {c.reference:.4g} got {c.predicted:.4g}"
              for c in run_experimental_validation() if not c.informational and not c.passed]
    assert not failed, "\n".join(failed)


def test_hall_refuses_tight_throat():
    # Hall diverges for Rc/Rt < 1.5 (axis M 1.30 at 0.625 vs 0.8 measured, Cuffel 1969)
    with pytest.raises(ValueError, match="diverges"):
        hall_line(1.4, 1.0, 0.625, axis_mach=1.02, n=10)
