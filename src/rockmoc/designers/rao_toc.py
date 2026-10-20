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
from rockmoc.designers.rao import design_rao_top
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

    # Build TOP at conical fraction first to get shape angles, then scale length
    # to length_fraction * L_ideal
    contour, meta = design_rao_top(
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach=Me,
        length_fraction=min(max(length_fraction, 0.6), 1.0),
        rc_over_rt=spec.throat_rc_over_rt,
        n_points=max(150, spec.n_characteristics * 2),
    )
    # Scale axial coordinate so total length = Lf * L_ideal
    L_target = length_fraction * L_ideal
    L_now = contour.length
    if L_now > 1e-12:
        scale = L_target / L_now
        from rockmoc.geometry.contour import Contour

        x0 = float(contour.x[0])
        x = x0 + scale * (contour.x - x0)
        contour = Contour(x=x, r=contour.r.copy(), theta=contour.theta, M=contour.M)

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
