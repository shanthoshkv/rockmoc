"""Post-solve dashboard: show nozzle + MOC + metrics after design."""

from __future__ import annotations

from typing import Any

from rockmoc.designers.base import DesignResult
from rockmoc.plot.moc_field import plot_moc_field
from rockmoc.plot.performance_plots import (
    plot_mach_pressure_along_wall,
    plot_nozzle_contour,
    plot_wall_profiles,
)


def display_design_result(
    result: DesignResult,
    *,
    show: bool = True,
    save_dir: str | None = None,
    prefix: str = "design",
    show_moc: bool = True,
    show_profiles: bool = True,
    show_thermo: bool = True,
) -> list[Any]:
    """Display (and optionally save) standard figures after a design solve.

    Saving is optional: pass ``save_dir`` to write PNGs; omit it to only display.

    Returns list of figure objects created.
    """
    figs: list[Any] = []
    metrics = {
        "kind": result.kind,
        "Me": f"{result.exit_mach_design:.3f}",
        "γ": f"{result.spec.gas.gamma:.3f}",
        "ε": f"{result.metrics.area_ratio:.3f}",
        "L/Rt": f"{result.metrics.length / result.spec.throat_radius:.3f}",
        "Cf,vac": f"{result.metrics.Cf_vac_1d:.4f}",
        "θ_max": f"{float(__import__('math').degrees(result.metrics.max_wall_angle)):.1f}°",
    }

    path = f"{save_dir}/{prefix}_contour.png" if save_dir else None
    figs.append(
        plot_nozzle_contour(
            result.contour,
            title=f"Nozzle contour — {result.kind.upper()}",
            metrics=metrics,
            show=show,
            save_path=path,
        )
    )

    if show_moc and result.net is not None:
        path = f"{save_dir}/{prefix}_moc_field.png" if save_dir else None
        figs.append(
            plot_moc_field(
                result.contour,
                result.net,
                title=f"MOC field (C+, C−) — {result.kind.upper()}",
                annotate=metrics,
                show=show,
                save_path=path,
            )
        )

    if show_profiles:
        path = f"{save_dir}/{prefix}_wall_profiles.png" if save_dir else None
        figs.append(
            plot_wall_profiles(
                result.contour,
                title=f"Wall profiles — {result.kind.upper()}",
                show=show,
                save_path=path,
            )
        )

    if show_thermo and result.contour.M is not None:
        path = f"{save_dir}/{prefix}_mach_pressure.png" if save_dir else None
        figs.append(
            plot_mach_pressure_along_wall(
                result.contour,
                result.spec.gas,
                title=f"Wall M & p/pc — {result.kind.upper()}",
                show=show,
                save_path=path,
            )
        )

    return figs


def display_production_result(
    prod: Any,
    *,
    show: bool = True,
    save_dir: str | None = None,
    prefix: str = "production",
) -> list[Any]:
    """Dashboard for ProductionResult: contour, MOC, BL, altitude."""
    from rockmoc.plot.performance_plots import (
        plot_bl_displacement,
        plot_performance_vs_altitude,
    )

    figs = display_design_result(
        prod.design, show=show, save_dir=save_dir, prefix=prefix, show_moc=True
    )
    if prod.bl is not None:
        path = f"{save_dir}/{prefix}_bl.png" if save_dir else None
        figs.append(
            plot_bl_displacement(
                prod.bl.inviscid,
                prod.bl.corrected,
                prod.bl.delta_star,
                show=show,
                save_path=path,
            )
        )
    if prod.altitude_sweep:
        path = f"{save_dir}/{prefix}_altitude.png" if save_dir else None
        figs.append(
            plot_performance_vs_altitude(
                prod.altitude_sweep, show=show, save_path=path
            )
        )
    return figs
