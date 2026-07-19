"""Structural typing for gases used in MOC (PerfectGas or VariableGammaGas)."""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from numpy.typing import ArrayLike, NDArray


@runtime_checkable
class MOCGas(Protocol):
    """Minimal gas API required by characteristic unit processes."""

    @property
    def gamma(self) -> float: ...

    @property
    def R(self) -> float: ...

    def prandtl_meyer(self, M: ArrayLike) -> NDArray | float: ...

    def mach_from_prandtl_meyer(self, nu: ArrayLike) -> NDArray | float: ...

    def mach_angle(self, M: ArrayLike) -> NDArray | float: ...

    def area_ratio(self, M: ArrayLike) -> NDArray | float: ...

    def mach_from_area_ratio(
        self, area_ratio: ArrayLike, *, branch: str = "supersonic"
    ) -> NDArray | float: ...

    def pe_pc_from_mach(self, M: ArrayLike) -> NDArray | float: ...

    def temperature_ratio(self, M: ArrayLike) -> NDArray | float: ...

    def pressure_ratio(self, M: ArrayLike) -> NDArray | float: ...

    def characteristic_velocity(self, T0: float) -> float: ...

    def thrust_coefficient(
        self,
        pe_pc: float,
        pa_pc: float,
        *,
        area_ratio: float | None = None,
        exit_mach: float | None = None,
    ) -> float: ...
