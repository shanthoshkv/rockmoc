"""Truncated ideal contour designer — high-fidelity path."""

from __future__ import annotations

from rockmoc.designers.axisym_hf import design_tic_hf
from rockmoc.designers.base import DesignResult, DesignSpec


def design_tic(spec: DesignSpec, *, store_net: bool = True) -> DesignResult:
    """Design a truncated ideal contour (fraction of ideal length)."""
    return design_tic_hf(spec, store_net=store_net)
