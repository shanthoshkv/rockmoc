"""Analysis: BL correction, throat IVL, atmosphere, comparisons, UQ. (lazy, PEP 562)."""

from __future__ import annotations

import importlib

_LAZY: dict[str, str] = {
    "BLCorrectionResult": "rockmoc.analysis.boundary_layer",
    "DesignComparison": "rockmoc.analysis.compare",
    "ThroatIVL": "rockmoc.analysis.throat",
    "altitude_from_pressure": "rockmoc.analysis.atmosphere",
    "circular_arc_throat_wall": "rockmoc.analysis.throat",
    "compare_nozzle_designs": "rockmoc.analysis.compare",
    "correct_wall_bl": "rockmoc.analysis.boundary_layer",
    "grid_convergence": "rockmoc.analysis.convergence",
    "mesh_certificate": "rockmoc.analysis.uq",
    "monte_carlo_mln": "rockmoc.analysis.uq",
    "pressure_at_altitude": "rockmoc.analysis.atmosphere",
    "richardson": "rockmoc.analysis.convergence",
    "sauer_throat_ivl": "rockmoc.analysis.throat",
    "sensitivity_one_at_a_time": "rockmoc.analysis.uq",
}

__all__ = [
    "BLCorrectionResult",
    "DesignComparison",
    "ThroatIVL",
    "altitude_from_pressure",
    "circular_arc_throat_wall",
    "compare_nozzle_designs",
    "correct_wall_bl",
    "grid_convergence",
    "mesh_certificate",
    "monte_carlo_mln",
    "pressure_at_altitude",
    "richardson",
    "sauer_throat_ivl",
    "sensitivity_one_at_a_time",
]

def __getattr__(name: str) -> object:
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError(f"module 'rockmoc.analysis' has no attribute {name!r}")
    value = getattr(importlib.import_module(module), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted({*globals(), *_LAZY})
