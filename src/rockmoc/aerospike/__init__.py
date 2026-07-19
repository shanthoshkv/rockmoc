"""Aerospike (plug) nozzle design — Lee / Vernacchia + base/cowl."""

from rockmoc.aerospike.base_cowl import (
    BasePressureModel,
    CowlDesign,
    design_primary_cowl,
    design_spike_with_cowl,
    estimate_base_pressure,
    truncated_spike_with_base,
)
from rockmoc.aerospike.lee_spike import (
    SpikeResult,
    design_spike,
    design_spike_truncated,
    majdalani_mach_from_er,
)

__all__ = [
    "BasePressureModel",
    "CowlDesign",
    "SpikeResult",
    "design_primary_cowl",
    "design_spike",
    "design_spike_truncated",
    "design_spike_with_cowl",
    "estimate_base_pressure",
    "majdalani_mach_from_er",
    "truncated_spike_with_base",
]
