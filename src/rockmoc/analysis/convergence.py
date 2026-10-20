"""Grid-convergence study with Richardson extrapolation and GCI (Roache 1994).

For three systematically refined meshes (N, rN, r²N) and a scalar f::

    p      = ln((f₁ − f₂) / (f₂ − f₃)) / ln r          observed order
    f_ext  = f₃ + (f₃ − f₂) / (r^p − 1)                Richardson estimate
    GCI₂₃  = F_s |(f₃ − f₂)/f₃| / (r^p − 1),  F_s = 1.25

Index 1 is the coarsest mesh.
"""

from __future__ import annotations

import math
from typing import Any, Literal

from rockmoc.moc.axisym_hf import design_axisym_ideal_hf, design_axisym_mln_hf
from rockmoc.thermo.gas_protocol import MOCGas

_METRICS = ("length_over_Rt", "area_ratio_geometric", "exit_wall_angle_deg", "max_wall_angle_deg")


def richardson(f: list[float], ratio: float = 2.0, fs: float = 1.25) -> dict[str, float | None]:
    """Observed order, extrapolated value and fine-grid GCI from three values (coarse → fine)."""
    f1, f2, f3 = f[-3:]
    d21, d32 = f2 - f1, f3 - f2
    tiny = 1e-9 * max(abs(f3), 1e-300)
    if abs(d32) <= tiny or abs(d21) <= tiny:  # converged to round-off
        return {"order": None, "extrapolated": f3, "gci_fine": 0.0}
    if d21 * d32 < 0.0:  # oscillatory convergence: no monotone order
        return {"order": None, "extrapolated": f3, "gci_fine": None}
    p = math.log(abs(d21 / d32)) / math.log(ratio)
    denom = ratio**p - 1.0
    return {
        "order": p,
        "extrapolated": f3 + d32 / denom,
        "gci_fine": abs(fs * (d32 / f3) / denom) if f3 != 0.0 else None,
    }


def grid_convergence(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    kind: Literal["mln", "ideal"] = "mln",
    n_base: int = 16,
    levels: int = 3,
    ratio: int = 2,
    delta: float = 1.0,
) -> dict[str, Any]:
    """Run ``levels`` refinements N = n_base·ratio^k and report order / extrapolation per metric."""
    if levels < 3:
        raise ValueError("levels must be ≥ 3 for an observed order")
    design = design_axisym_mln_hf if kind == "mln" else design_axisym_ideal_hf
    rows: list[dict[str, float]] = []
    for k in range(levels):
        n = n_base * ratio**k
        res = design(gas, throat_radius=throat_radius, exit_mach=exit_mach, n_chars=n,
                     delta=delta, match_exit_area=False, store_net=True)
        w = res.net.wall if res.net is not None else []
        last = math.degrees(math.atan2(w[-1].r - w[-2].r, w[-1].x - w[-2].x)) if len(w) > 1 else 0.0
        s = res.summary()
        rows.append({
            "n_chars": n,
            "length_over_Rt": s["length_over_Rt"],
            "area_ratio_geometric": s["area_ratio_geometric"],
            "exit_wall_angle_deg": abs(last),
            "max_wall_angle_deg": s["max_wall_angle_deg"],
        })
    stats = {m: richardson([r[m] for r in rows], float(ratio)) for m in _METRICS}
    return {"kind": kind, "exit_mach": exit_mach, "rows": rows, "richardson": stats}
