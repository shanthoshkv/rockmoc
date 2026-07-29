"""
rockmoc — Industry-oriented MOC rocket nozzle design library.

Production pipeline::

    from rockmoc import DesignSpec, PerfectGas, design_production
    r = design_production(DesignSpec(PerfectGas(1.25), throat_radius=0.05, exit_mach=3.0))
    print(r.summary())
"""

from rockmoc._version import __version__
from rockmoc.aerospike import (
    SpikeResult,
    design_primary_cowl,
    design_spike,
    design_spike_truncated,
    design_spike_with_cowl,
    truncated_spike_with_base,
)
from rockmoc.analysis import (
    compare_nozzle_designs,
    correct_wall_bl,
    mesh_certificate,
    monte_carlo_mln,
    pressure_at_altitude,
    sauer_throat_ivl,
    sensitivity_one_at_a_time,
)
from rockmoc.designers import (
    DesignResult,
    DesignSpec,
    ProductionResult,
    design_conical,
    design_ideal,
    design_mln,
    design_production,
    design_rao,
    design_rao_control,
    design_rao_toc,
    design_tic,
)
from rockmoc.student import (
    KNOWN_LIMITS,
    comparison_to_csv,
    comparison_to_latex,
    comparison_to_markdown,
    describe_example_gases,
    example_gas,
    explainer,
    format_design_report,
    generate_lab_report,
    list_example_gases,
)
from rockmoc.exceptions import (
    ConvergenceError,
    DesignError,
    MeshError,
    RockmocError,
    ThermoError,
)
from rockmoc.geometry import Contour
from rockmoc.moc.axisym_hf import mesh_convergence_axisym
from rockmoc.moc.planar_nasa import PlanarNasaResult, design_planar_mln_nasa, mesh_convergence_planar
from rockmoc.performance import (
    field_integrated_cf,
    integrated_performance,
    off_design_altitude_sweep,
)
from rockmoc.thermo import (
    PerfectGas,
    VariableGammaGas,
    build_multi_gamma_table,
    build_perfect_gas_table,
)
from rockmoc.validate import run_all_validation, validation_report_markdown

__all__ = [
    "__version__",
    "Contour",
    "DesignResult",
    "DesignSpec",
    "PerfectGas",
    "PlanarNasaResult",
    "ProductionResult",
    "SpikeResult",
    "VariableGammaGas",
    "RockmocError",
    "DesignError",
    "ThermoError",
    "MeshError",
    "ConvergenceError",
    "design_mln",
    "design_ideal",
    "design_tic",
    "design_rao",
    "design_rao_control",
    "design_rao_toc",
    "design_conical",
    "design_production",
    "design_planar_mln_nasa",
    "mesh_convergence_planar",
    "mesh_convergence_axisym",
    "design_spike",
    "design_spike_truncated",
    "design_spike_with_cowl",
    "design_primary_cowl",
    "truncated_spike_with_base",
    "correct_wall_bl",
    "sauer_throat_ivl",
    "compare_nozzle_designs",
    "comparison_to_markdown",
    "comparison_to_latex",
    "comparison_to_csv",
    "example_gas",
    "list_example_gases",
    "describe_example_gases",
    "explainer",
    "format_design_report",
    "generate_lab_report",
    "KNOWN_LIMITS",
    "integrated_performance",
    "field_integrated_cf",
    "off_design_altitude_sweep",
    "pressure_at_altitude",
    "monte_carlo_mln",
    "mesh_certificate",
    "sensitivity_one_at_a_time",
    "build_perfect_gas_table",
    "build_multi_gamma_table",
    "run_all_validation",
    "validation_report_markdown",
]
