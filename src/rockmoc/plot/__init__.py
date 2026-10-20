"""Dark-themed plotting for rockmoc (optional matplotlib)."""

from rockmoc.plot.analysis_plots import plot_grid_convergence
from rockmoc.plot.dashboard import display_design_result, display_production_result
from rockmoc.plot.moc_field import plot_moc_field
from rockmoc.plot.performance_plots import (
    plot_aerospike,
    plot_area_ratio_sweep,
    plot_bl_displacement,
    plot_contour_overlay,
    plot_design_comparison,
    plot_gamma_sweep,
    plot_mach_pressure_along_wall,
    plot_nozzle_contour,
    plot_of_ratio_sweep,
    plot_performance_vs_altitude,
    plot_wall_profiles,
    plot_wall_properties,
)
from rockmoc.plot.sweeps import (
    altitude_sweep_rows,
    area_ratio_sweep_rows,
    gamma_sweep_rows,
    of_ratio_sweep_rows,
)
from rockmoc.plot.theme import COLORS, apply_dark_theme, themed


def plot_moc_field_interactive(*args, **kwargs):  # noqa: ANN002, ANN003
    """Plotly MOC field (lazy import; needs ``rockmoc[viz]``)."""
    from rockmoc.plot.interactive import plot_moc_field_interactive as _f

    return _f(*args, **kwargs)


__all__ = [
    "COLORS",
    "apply_dark_theme",
    "altitude_sweep_rows",
    "area_ratio_sweep_rows",
    "display_design_result",
    "display_production_result",
    "gamma_sweep_rows",
    "of_ratio_sweep_rows",
    "plot_aerospike",
    "plot_area_ratio_sweep",
    "plot_bl_displacement",
    "plot_contour_overlay",
    "plot_design_comparison",
    "plot_grid_convergence",
    "plot_gamma_sweep",
    "plot_mach_pressure_along_wall",
    "plot_moc_field",
    "plot_moc_field_interactive",
    "plot_nozzle_contour",
    "plot_of_ratio_sweep",
    "plot_performance_vs_altitude",
    "plot_wall_profiles",
    "plot_wall_properties",
    "themed",
]
