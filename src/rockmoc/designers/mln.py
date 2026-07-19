"""Minimum-length nozzle (MLN) designer — high-fidelity axisymmetric MOC."""

from __future__ import annotations

from rockmoc.designers.axisym_hf import design_mln_hf
from rockmoc.designers.base import DesignResult, DesignSpec


def design_mln(spec: DesignSpec, *, store_net: bool = True) -> DesignResult:
    """Design a minimum-length rocket nozzle contour (axisymmetric MOC).

    Uses high-fidelity unit-process MOC (Zucrow/Hoffman-class). Corner angle
    starts from planar theory θ_max = ν_e/2 and is lightly iterated for
    axisymmetric flow. Geometric exit area ratio is an **output** (not forced).

    For the exact NASA **planar** TM-X-1502 algorithm, use
    :func:`rockmoc.design_planar_mln_nasa` instead.
    """
    return design_mln_hf(spec, store_net=store_net)
