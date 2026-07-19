"""Fair comparison of MLN / ideal / TIC / Rao on the same design point."""

from __future__ import annotations

from dataclasses import dataclass

from rockmoc.designers import DesignSpec, design_ideal, design_mln, design_tic
from rockmoc.designers.rao import design_rao
from rockmoc.performance.integrated import integrated_performance


@dataclass
class DesignComparison:
    rows: list[dict]

    def summary_table(self) -> str:
        hdr = f"{'kind':10} {'L/Rt':>8} {'eps':>8} {'Cf_vac':>8} {'eta_n':>8} {'th_max':>8}"
        lines = [hdr, "-" * len(hdr)]
        for r in self.rows:
            lines.append(
                f"{r['kind']:10} {r['L_over_Rt']:8.3f} {r['area_ratio']:8.3f} "
                f"{r['Cf_vac']:8.4f} {r['eta_n']:8.4f} {r['theta_max_deg']:8.2f}"
            )
        return "\n".join(lines)


def compare_nozzle_designs(
    spec: DesignSpec,
    *,
    rao_length_fraction: float = 0.80,
    tic_fraction: float = 0.80,
    pc: float = 5e6,
    Tc: float = 3000.0,
) -> DesignComparison:
    """Design MLN, ideal, TIC, Rao at the same Me/γ/Rt and compare metrics."""
    spec.validate()
    Me = spec.resolved_exit_mach()
    Rt = spec.throat_radius
    gas = spec.gas

    designs = {
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
                "exit_radius_m": res.metrics.exit_radius,
            }
        )
    return DesignComparison(rows=rows)
