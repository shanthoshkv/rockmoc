"""rockmoc — Industry-oriented MOC rocket nozzle design library.

Production pipeline::

    from rockmoc import DesignSpec, PerfectGas, design_production
    r = design_production(DesignSpec(PerfectGas(1.25), throat_radius=0.05, exit_mach=3.0))
    print(r.summary())

Public names are imported lazily on first access (PEP 562), so
``import rockmoc`` stays fast; ``from rockmoc import X`` works as before.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

from rockmoc._version import __version__

if TYPE_CHECKING:  # static analysers / IDEs see the real imports
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
    from rockmoc.exceptions import (
        ConvergenceError,
        DesignError,
        MeshError,
        RockmocError,
        ThermoError,
    )
    from rockmoc.geometry import Contour
    from rockmoc.moc.axisym_hf import mesh_convergence_axisym
    from rockmoc.moc.planar_nasa import (
        PlanarNasaResult,
        design_planar_mln_nasa,
        mesh_convergence_planar,
    )
    from rockmoc.performance import (
        field_integrated_cf,
        integrated_performance,
        off_design_altitude_sweep,
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
    from rockmoc.thermo import (
        PerfectGas,
        VariableGammaGas,
        build_multi_gamma_table,
        build_perfect_gas_table,
    )
    from rockmoc.validate import run_all_validation, validation_report_markdown

_LAZY: dict[str, str] = {
    "Contour": "rockmoc.geometry",
    "ConvergenceError": "rockmoc.exceptions",
    "DesignError": "rockmoc.exceptions",
    "DesignResult": "rockmoc.designers",
    "DesignSpec": "rockmoc.designers",
    "KNOWN_LIMITS": "rockmoc.student",
    "MeshError": "rockmoc.exceptions",
    "PerfectGas": "rockmoc.thermo",
    "PlanarNasaResult": "rockmoc.moc.planar_nasa",
    "ProductionResult": "rockmoc.designers",
    "RockmocError": "rockmoc.exceptions",
    "SpikeResult": "rockmoc.aerospike",
    "ThermoError": "rockmoc.exceptions",
    "VariableGammaGas": "rockmoc.thermo",
    "build_multi_gamma_table": "rockmoc.thermo",
    "build_perfect_gas_table": "rockmoc.thermo",
    "compare_nozzle_designs": "rockmoc.analysis",
    "comparison_to_csv": "rockmoc.student",
    "comparison_to_latex": "rockmoc.student",
    "comparison_to_markdown": "rockmoc.student",
    "correct_wall_bl": "rockmoc.analysis",
    "describe_example_gases": "rockmoc.student",
    "design_conical": "rockmoc.designers",
    "design_ideal": "rockmoc.designers",
    "design_mln": "rockmoc.designers",
    "design_planar_mln_nasa": "rockmoc.moc.planar_nasa",
    "design_primary_cowl": "rockmoc.aerospike",
    "design_production": "rockmoc.designers",
    "design_rao": "rockmoc.designers",
    "design_rao_control": "rockmoc.designers",
    "design_rao_toc": "rockmoc.designers",
    "design_spike": "rockmoc.aerospike",
    "design_spike_truncated": "rockmoc.aerospike",
    "design_spike_with_cowl": "rockmoc.aerospike",
    "design_tic": "rockmoc.designers",
    "example_gas": "rockmoc.student",
    "explainer": "rockmoc.student",
    "field_integrated_cf": "rockmoc.performance",
    "format_design_report": "rockmoc.student",
    "generate_lab_report": "rockmoc.student",
    "integrated_performance": "rockmoc.performance",
    "list_example_gases": "rockmoc.student",
    "mesh_certificate": "rockmoc.analysis",
    "mesh_convergence_axisym": "rockmoc.moc.axisym_hf",
    "mesh_convergence_planar": "rockmoc.moc.planar_nasa",
    "monte_carlo_mln": "rockmoc.analysis",
    "off_design_altitude_sweep": "rockmoc.performance",
    "pressure_at_altitude": "rockmoc.analysis",
    "run_all_validation": "rockmoc.validate",
    "sauer_throat_ivl": "rockmoc.analysis",
    "sensitivity_one_at_a_time": "rockmoc.analysis",
    "truncated_spike_with_base": "rockmoc.aerospike",
    "validation_report_markdown": "rockmoc.validate",
}

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

def __getattr__(name: str) -> object:
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError(f"module 'rockmoc' has no attribute {name!r}")
    value = getattr(importlib.import_module(module), name)
    globals()[name] = value  # cache: later lookups bypass __getattr__
    return value


def __dir__() -> list[str]:
    return sorted({*globals(), *_LAZY})
