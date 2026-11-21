"""Public designers wrapping high-fidelity axisymmetric MOC."""

from __future__ import annotations

import math

import numpy as np

from rockmoc.designers.base import DesignResult, DesignSpec
from rockmoc.geometry.contour import Contour
from rockmoc.moc.axisym_hf import (
    design_axisym_ideal_hf,
    design_axisym_mln_hf,
    design_axisym_tic_hf,
)
from rockmoc.performance.metrics import compute_metrics


def _densify_throat_arc(c: Contour, theta_attach_deg: float, n_arc: int = 40) -> Contour:
    """Replace the few MOC wall nodes on the circular throat arc by ``n_arc`` exact arc points.

    The net hits the arc only where its C+ lines land (4–5 nodes for n≈40), so a
    polyline/CAD fit of the contour carries 3–6° slope kinks right at the throat;
    CFD showed the resulting waves as a 1.5–3 % Mach error over the first 2 Rt.
    """
    th, x, r = np.asarray(c.theta), np.asarray(c.x), np.asarray(c.r)
    above = np.nonzero(th > math.radians(theta_attach_deg) + 1e-9)[0]  # first node past the arc
    k = int(above[0]) - 1 if above.size else 0
    if k < 1 or th[k] <= 0.0:
        return c
    Rt, x0 = float(r[0]), float(x[0])
    Rc = (x[k] - x0) / math.sin(th[k])  # effective arc radius (after exit-area matching)
    phi = np.linspace(0.0, th[k], n_arc)
    xa, ra = x0 + Rc * np.sin(phi), Rt + Rc * (1.0 - np.cos(phi))
    ra[-1] = r[k]
    M = (
        None
        if c.M is None
        else np.concatenate([np.interp(xa[:-1], x[: k + 1], np.asarray(c.M)[: k + 1]), c.M[k:]])
    )
    return Contour(
        x=np.concatenate([xa[:-1], x[k:]]),
        r=np.concatenate([ra[:-1], r[k:]]),
        theta=np.concatenate([phi[:-1], th[k:]]),
        M=M,
    )


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
    contour = hf.contour  # sharp-corner throat: no arc to densify
    metrics = compute_metrics(
        contour,
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach_design=Me,
    )
    return DesignResult(
        kind="mln",
        spec=spec,
        exit_mach_design=Me,
        contour=contour,
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
    contour = _densify_throat_arc(hf.contour, hf.summary()["theta_attach_deg"])
    metrics = compute_metrics(
        contour,
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach_design=Me,
    )
    return DesignResult(
        kind="ideal",
        spec=spec,
        exit_mach_design=Me,
        contour=contour,
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
    contour = _densify_throat_arc(hf.contour, hf.summary()["theta_attach_deg"])
    metrics = compute_metrics(
        contour,
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach_design=Me,
    )
    return DesignResult(
        kind="tic",
        spec=spec,
        exit_mach_design=Me,
        contour=contour,
        metrics=metrics,
        net=hf.net,
        meta={**hf.summary(), "designer": "axisym_hf"},
    )
