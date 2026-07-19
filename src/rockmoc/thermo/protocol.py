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

from typing import Protocol, runtime_checkable

from numpy.typing import ArrayLike, NDArray


@runtime_checkable
class GasModel(Protocol):
    """Minimal thermo interface required by MOC designers."""

    @property
    def gamma(self) -> float:
        """Effective ratio of specific heats for PM / area-Mach (design γ)."""
        ...

    @property
    def R(self) -> float:
        """Specific gas constant [J/(kg·K)]."""
        ...

    def prandtl_meyer(self, M: ArrayLike) -> NDArray | float:
        """Prandtl–Meyer function ν(M) [radians], M ≥ 1."""
        ...

    def mach_from_prandtl_meyer(self, nu: ArrayLike) -> NDArray | float:
        """Inverse Prandtl–Meyer M(ν)."""
        ...

    def mach_angle(self, M: ArrayLike) -> NDArray | float:
        """Mach angle μ = arcsin(1/M) [radians]."""
        ...

    def area_ratio(self, M: ArrayLike) -> NDArray | float:
        """Isentropic A/A* for Mach M."""
        ...

    def mach_from_area_ratio(
        self, area_ratio: ArrayLike, *, branch: str = "supersonic"
    ) -> NDArray | float:
        """Invert A/A* → M."""
        ...

    def pressure_ratio(self, M: ArrayLike) -> NDArray | float:
        """Stagnation / static pressure p0/p."""
        ...

    def temperature_ratio(self, M: ArrayLike) -> NDArray | float:
        """Stagnation / static temperature T0/T."""
        ...
