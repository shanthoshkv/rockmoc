"""Public designers wrapping high-fidelity axisymmetric MOC."""

from __future__ import annotations

from rockmoc.designers.base import DesignResult, DesignSpec
from rockmoc.moc.axisym_hf import (
    design_axisym_ideal_hf,
    design_axisym_mln_hf,
    design_axisym_tic_hf,
)
from rockmoc.performance.metrics import compute_metrics


def design_mln_hf(spec: DesignSpec, *, store_net: bool = True) -> DesignResult:
    """High-fidelity axisymmetric/planar MLN (unit-process MOC)."""
    spec.validate()
    Me = spec.resolved_exit_mach()
    # Respect user n_characteristics (production: prefer ≥ 80; tests may use less)
    n = max(8, int(spec.n_characteristics))

    hf = design_axisym_mln_hf(
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach=Me,
        n_chars=n,
        throat_mach=spec.throat_mach,
        delta=spec.delta,
        match_exit_area=True,
        store_net=store_net,
    )
    metrics = compute_metrics(
        hf.contour,
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach_design=Me,
    )
    return DesignResult(
        kind="mln",
        spec=spec,
        exit_mach_design=Me,
        contour=hf.contour,
        metrics=metrics,
        net=hf.net,
        meta={**hf.summary(), "designer": "axisym_hf"},
    )


def design_ideal_hf(spec: DesignSpec, *, store_net: bool = True) -> DesignResult:
    """High-fidelity ideal-like axisymmetric contour."""
    spec.validate()
    Me = spec.resolved_exit_mach()
    n = max(8, int(spec.n_characteristics))

    hf = design_axisym_ideal_hf(
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach=Me,
        n_chars=n,
        throat_rc_over_rt=spec.throat_rc_over_rt,
        throat_mach=spec.throat_mach,
        delta=spec.delta,
        match_exit_area=True,
        store_net=store_net,
    )
    metrics = compute_metrics(
        hf.contour,
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach_design=Me,
    )
    return DesignResult(
        kind="ideal",
        spec=spec,
        exit_mach_design=Me,
        contour=hf.contour,
        metrics=metrics,
        net=hf.net,
        meta={**hf.summary(), "designer": "axisym_hf"},
    )


def design_tic_hf(spec: DesignSpec, *, store_net: bool = True) -> DesignResult:
    """High-fidelity truncated ideal contour."""
    spec.validate()
    Me = spec.resolved_exit_mach()
    frac = spec.truncate_fraction if spec.truncate_fraction is not None else 0.80
    n = max(8, int(spec.n_characteristics))

    hf = design_axisym_tic_hf(
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach=Me,
        truncate_fraction=frac,
        n_chars=n,
        throat_rc_over_rt=spec.throat_rc_over_rt,
        throat_mach=spec.throat_mach,
        delta=spec.delta,
        match_exit_area=True,
        store_net=store_net,
    )
    metrics = compute_metrics(
        hf.contour,
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach_design=Me,
    )
    return DesignResult(
        kind="tic",
        spec=spec,
        exit_mach_design=Me,
        contour=hf.contour,
        metrics=metrics,
        net=hf.net,
        meta={**hf.summary(), "designer": "axisym_hf"},
    )
