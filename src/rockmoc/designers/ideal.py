"""Ideal contour designer — high-fidelity axisymmetric MOC."""

from __future__ import annotations

from rockmoc.designers.axisym_hf import design_ideal_hf
from rockmoc.designers.base import DesignResult, DesignSpec


def design_ideal(spec: DesignSpec, *, store_net: bool = True) -> DesignResult:
    """Design an ideal-like axisymmetric contour (circular throat + elongated wall)."""
    return design_ideal_hf(spec, store_net=store_net)
