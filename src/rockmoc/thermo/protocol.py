"""Gas model protocol for constant-γ and chemistry-backed properties.

Design philosophy
-----------------
MOC unit processes need local thermodynamic response along isentropic
paths. Industry practice usually freezes composition at the throat and
uses a single effective γ (RocketCEA "frozen"). Higher fidelity uses
isentropic equilibrium tables (variable γ_s).

Backends implement this protocol; the core MOC never imports CEA/Cantera
directly.
"""

from __future__ import annotations

from rockmoc.thermo.gas_protocol import MOCGas

# Single gas interface for the whole library (kept under its historic name).
GasModel = MOCGas

__all__ = ["GasModel"]
