"""Thrust-optimized contour (TOC) closer to Rao's method.

Uses the MOC ideal kernel length as the "100%" reference, then designs a
TOP-style bell at a prescribed length fraction with Rao θ_n, θ_e, and
matches the design area ratio. Not a full calculus-of-variations solver,
but uses MOC ideal length rather than only a 15° cone reference.

No CFD.
"""

from __future__ import annotations

from rockmoc.designers.base import DesignResult, DesignSpec
from rockmoc.designers.ideal import design_ideal
from rockmoc.designers.rao import conical_length, design_rao_top
from rockmoc.performance.metrics import compute_metrics


def design_rao_toc(
    spec: DesignSpec,
    *,
    length_fraction: float = 0.80,
    store_net: bool = False,
) -> DesignResult:
    """Rao TOC: length = fraction of *MOC ideal* length; TOP parabola wall.

    Parameters
    ----------
    length_fraction:
        Fraction of the ideal MOC nozzle length (not conical). Typical 0.6–0.9.
    """
    spec.validate()
    Me = spec.resolved_exit_mach()
    # Ideal MOC length as reference
    ideal = design_ideal(spec, store_net=False)
    L_ideal = ideal.metrics.length
    if L_ideal <= 0:
        L_ideal = 10.0 * spec.throat_radius

    # TOP bell built directly at L = length_fraction · L_ideal: circular throat arc kept,
    # Rao angles taken at the equivalent 15°-cone length fraction (clipped to the
    # 0.6–1.0 chart range). v0.9 stretched a TOP contour in x afterwards, which
    # left θ unscaled (stated θ_e 10.7° vs 5.6° real) and turned the Rc = 1.5 Rt
    # throat arc into an ellipse with Rc ≈ 5.6 Rt.
    L_target = length_fraction * L_ideal
    eps = float(spec.gas.area_ratio(Me))
    Lf_cone = L_target / conical_length(spec.throat_radius, eps, 15.0)
    contour, meta = design_rao_top(
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach=Me,
        length_fraction=min(max(Lf_cone, 0.6), 1.0),
        rc_over_rt=spec.throat_rc_over_rt,
        n_points=max(150, spec.n_characteristics * 2),
        length_m=L_target,
    )
    meta["length_fraction_of_cone"] = Lf_cone

    metrics = compute_metrics(
        contour,
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach_design=Me,
    )
    meta = {
        **meta,
        "method": "rao_toc_moc_ideal_length",
        "L_ideal_moc_m": L_ideal,
        "length_fraction_of_ideal": length_fraction,
        "length_m": metrics.length,
        "designer": "rao_toc",
        "reference": "Rao TOP + MOC ideal length reference (no CFD)",
    }
    return DesignResult(
        kind="rao",
        spec=spec,
        exit_mach_design=Me,
        contour=contour,
        metrics=metrics,
        net=None,
        meta=meta,
    )
