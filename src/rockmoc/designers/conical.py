"""Conical nozzle baseline — textbook comparison shape."""

from __future__ import annotations

import math

import numpy as np

from rockmoc.designers.base import DesignResult, DesignSpec
from rockmoc.exceptions import DesignError
from rockmoc.geometry.contour import Contour
from rockmoc.performance.metrics import compute_metrics


def design_conical(
    spec: DesignSpec,
    *,
    half_angle_deg: float = 15.0,
    n_points: int = 80,
    store_net: bool = False,
) -> DesignResult:
    """Axisymmetric conical diverging nozzle (constant wall half-angle).

    The exit radius is set from 1-D isentropic area ratio at the design Me.
    Length follows pure geometry::

        L = (R_e − R_t) / tan(θ_c)

    This is the standard homework baseline used to judge MLN / TIC / Rao gains.

    Parameters
    ----------
    half_angle_deg:
        Constant wall half-angle [degrees]. Classic teaching value is **15°**.
    n_points:
        Number of wall samples (≥ 3).
    store_net:
        Unused (cones have no MOC net); accepted for API symmetry.
    """
    _ = store_net
    spec.validate()
    if half_angle_deg <= 0.0 or half_angle_deg >= 45.0:
        raise DesignError(
            f"half_angle_deg must be in (0, 45); got {half_angle_deg}",
            tip="Use a textbook cone angle such as 12–18° (default 15°).",
            example="design_conical(spec, half_angle_deg=15.0)",
        )
    if n_points < 3:
        raise DesignError(
            "n_points must be ≥ 3",
            tip="Increase wall sampling for a smooth CSV export.",
            example="design_conical(spec, n_points=80)",
        )

    Rt = float(spec.throat_radius)
    Me = float(spec.resolved_exit_mach())
    gas = spec.gas
    eps = float(gas.area_ratio(Me))
    Re = Rt * math.sqrt(eps)
    th = math.radians(float(half_angle_deg))
    L = (Re - Rt) / math.tan(th)
    if L <= 0.0:
        raise DesignError(
            "Cone length non-positive — check Me / area_ratio / throat_radius",
            tip="Design Me must be > 1 so exit radius exceeds the throat.",
        )

    x = np.linspace(0.0, L, n_points)
    r = Rt + x * math.tan(th)
    # Force exact end radius
    r[-1] = Re
    theta = np.full(n_points, th)
    # Edge Mach from local area ratio (1-D isentropic along cone)
    M = np.empty(n_points)
    for i in range(n_points):
        ar = max((r[i] / Rt) ** 2, 1.0 + 1e-12)
        try:
            M[i] = float(gas.mach_from_area_ratio(ar, branch="supersonic"))
        except Exception:
            M[i] = 1.01 + (Me - 1.01) * (i / max(n_points - 1, 1))
    M[0] = max(float(spec.throat_mach), 1.01)
    M[-1] = Me

    contour = Contour(x=x, r=r, theta=theta, M=M)
    metrics = compute_metrics(
        contour, gas, throat_radius=Rt, exit_mach_design=Me
    )
    return DesignResult(
        kind="conical",  # type: ignore[arg-type]
        spec=spec,
        exit_mach_design=Me,
        contour=contour,
        metrics=metrics,
        net=None,
        meta={
            "method": "conical_constant_half_angle",
            "half_angle_deg": float(half_angle_deg),
            "half_angle_rad": th,
            "area_ratio_1d": eps,
            "explainer": "conical",
            "limitations": (
                "Straight-wall cone; large exit angle → divergence loss. "
                "Teaching baseline only."
            ),
        },
    )
