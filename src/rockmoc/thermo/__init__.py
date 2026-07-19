"""Thermodynamic relations for perfect-gas and chemistry-backed nozzle design."""

from rockmoc.thermo.perfect_gas import (
    PerfectGas,
    area_ratio,
    mach_angle,
    mach_from_area_ratio,
    mach_from_prandtl_meyer,
    prandtl_meyer,
)
from rockmoc.thermo.expansion_table import (
    ExpansionTable,
    build_multi_gamma_table,
    build_perfect_gas_table,
)
from rockmoc.thermo.protocol import GasModel
from rockmoc.thermo.variable_gamma import VariableGammaGas

# Optional backends imported lazily by users:
#   from rockmoc.thermo.rocketcea_backend import RocketCEAGas
#   from rockmoc.thermo.cantera_backend import CanteraGas
#   from rockmoc.thermo.coolprop_backend import CoolPropGas
#   from rockmoc.thermo.expansion_table import build_cea_frozen_table

__all__ = [
    "ExpansionTable",
    "GasModel",
    "PerfectGas",
    "VariableGammaGas",
    "area_ratio",
    "build_multi_gamma_table",
    "build_perfect_gas_table",
    "mach_angle",
    "mach_from_area_ratio",
    "mach_from_prandtl_meyer",
    "prandtl_meyer",
]
