"""Fair comparison of MLN / ideal / TIC / Rao / conical on the same design point."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rockmoc.designers import DesignSpec, design_ideal, design_mln, design_tic
from rockmoc.designers.conical import design_conical
from rockmoc.designers.rao import design_rao
from rockmoc.performance.integrated import integrated_performance


@dataclass
class DesignComparison:
    """Table of design metrics for side-by-side export."""

    rows: list[dict]

    def summary_table(self) -> str:
        hdr = (
            f"{'kind':10} {'L/Rt':>8} {'eps':>8} {'Cf_vac':>8} "
            f"{'eta_n':>8} {'th_max':>8} {'th_e':>8}"
        )
        lines = [hdr, "-" * len(hdr)]
        for r in self.rows:
            lines.append(
                f"{r['kind']:10} {r['L_over_Rt']:8.3f} {r['area_ratio']:8.3f} "
                f"{r['Cf_vac']:8.4f} {r['eta_n']:8.4f} "
                f"{r['theta_max_deg']:8.2f} {r.get('theta_exit_deg', 0.0):8.2f}"
            )
        return "\n".join(lines)

    def to_markdown(self, **kwargs: Any) -> str:
        from rockmoc.student.tables import comparison_to_markdown

        return comparison_to_markdown(self, **kwargs)

    def to_latex(self, **kwargs: Any) -> str:
        from rockmoc.student.tables import comparison_to_latex

        return comparison_to_latex(self, **kwargs)

    def to_csv(self, path: str | Path | None = None, **kwargs: Any) -> str:
        from rockmoc.student.tables import comparison_to_csv

        return comparison_to_csv(self, path, **kwargs)

    def to_string(self) -> str:
        return self.summary_table()


def compare_nozzle_designs(
    spec: DesignSpec,
    *,
    rao_length_fraction: float = 0.80,
    tic_fraction: float = 0.80,
    cone_half_angle_deg: float = 15.0,
    include_conical: bool = True,
    pc: float = 5e6,
    Tc: float = 3000.0,
) -> DesignComparison:
    """Design MLN, ideal, TIC, Rao, and optional conical baseline; compare metrics."""
    spec.validate()
    Me = spec.resolved_exit_mach()
    Rt = spec.throat_radius
    gas = spec.gas

    designs: dict[str, Any] = {
        "mln": design_mln(spec, store_net=False),
        "ideal": design_ideal(spec, store_net=False),
        "tic": design_tic(
            DesignSpec(
                gas=gas,
                throat_radius=Rt,
                exit_mach=Me,
                n_characteristics=spec.n_characteristics,
                truncate_fraction=tic_fraction,
                throat_rc_over_rt=spec.throat_rc_over_rt,
            ),
            store_net=False,
        ),
        "rao": design_rao(spec, length_fraction=rao_length_fraction),
    }
    if include_conical:
        designs["conical"] = design_conical(
            spec, half_angle_deg=cone_half_angle_deg, n_points=60
        )

    rows = []
    for kind, res in designs.items():
        perf = integrated_performance(
            res.contour,
            gas,
            throat_radius=Rt,
            exit_mach=Me,
            pc=pc,
            Tc=Tc,
            pa=0.0,
        )
        th_max = res.metrics.max_wall_angle
        th_e = (
            float(abs(res.contour.theta[-1]))
            if res.contour.theta is not None
            else 0.0
        )
        rows.append(
            {
                "kind": kind,
                "L_m": res.metrics.length,
                "L_over_Rt": res.metrics.length / Rt,
                "area_ratio": res.metrics.area_ratio,
                "Cf_vac": perf.Cf_vac,
                "Isp_vac_s": perf.Isp_vac,
                "eta_n": perf.nozzle_efficiency,
                "eta_div": perf.divergence_efficiency,
                "theta_max_deg": float(__import__("math").degrees(th_max)),
                "theta_exit_deg": float(__import__("math").degrees(th_e)),
                "exit_radius_m": res.metrics.exit_radius,
            }
        )
    return DesignComparison(rows=rows)
