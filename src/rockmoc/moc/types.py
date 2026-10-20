"""Data structures for the characteristic mesh."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum, auto

import numpy as np
from numpy.typing import NDArray

from rockmoc.exceptions import ThermoError
from rockmoc.thermo.gas_protocol import MOCGas


class PointKind(Enum):
    """Classification of a mesh node."""

    INTERIOR = auto()
    AXIS = auto()
    WALL = auto()
    FAN = auto()
    THROAT = auto()
    EXIT = auto()


@dataclass(slots=True, frozen=True)
class FlowPoint:
    """State and position of a single characteristic-mesh node.

    Angles ``theta``, ``nu``, and ``mu`` are stored in **radians**.
    Coordinates ``x``, ``r`` are in meters (SI).
    ``gamma_local`` is optional local γ for variable-γ tables (else NaN).
    """

    x: float
    r: float
    theta: float
    nu: float
    M: float
    mu: float
    kind: PointKind = PointKind.INTERIOR
    id: int = -1
    gamma_local: float = float("nan")

    def copy(self, **kwargs: object) -> FlowPoint:
        return replace(self, **kwargs)  # type: ignore[arg-type]

    @classmethod
    def from_mach_theta(
        cls,
        x: float,
        r: float,
        M: float,
        theta: float,
        gas: MOCGas,
        *,
        kind: PointKind = PointKind.INTERIOR,
        id: int = -1,
    ) -> FlowPoint:
        if M < 1.0:
            raise ThermoError(f"FlowPoint requires M ≥ 1; got {M}")
        nu = float(gas.prandtl_meyer(M))
        mu = float(gas.mach_angle(M))
        g_loc = float(gas.gamma_at_mach(M))
        return cls(
            x=x, r=r, theta=theta, nu=nu, M=M, mu=mu, kind=kind, id=id, gamma_local=g_loc
        )

    @classmethod
    def from_nu_theta(
        cls,
        x: float,
        r: float,
        nu: float,
        theta: float,
        gas: MOCGas,
        *,
        kind: PointKind = PointKind.INTERIOR,
        id: int = -1,
    ) -> FlowPoint:
        M = float(gas.mach_from_prandtl_meyer(nu))
        mu = float(gas.mach_angle(M))
        g_loc = float(gas.gamma_at_mach(M))
        return cls(
            x=x, r=r, theta=theta, nu=nu, M=M, mu=mu, kind=kind, id=id, gamma_local=g_loc
        )


@dataclass
class CharNet:
    """Collection of characteristic points and connectivity."""

    points: list[FlowPoint] = field(default_factory=list)
    # Optional index grids for structured nets (row = characteristic index)
    interior: list[list[FlowPoint]] = field(default_factory=list)
    axis: list[FlowPoint] = field(default_factory=list)
    wall: list[FlowPoint] = field(default_factory=list)
    fan: list[FlowPoint] = field(default_factory=list)
    # C− characteristic polylines (upstream → downstream); ``interior`` rows are C+
    c_minus: list[list[FlowPoint]] = field(default_factory=list)

    def add(self, p: FlowPoint) -> FlowPoint:
        object.__setattr__(p, "id", len(self.points))  # bookkeeping only; FlowPoint is frozen
        self.points.append(p)
        return p

    def as_arrays(self) -> dict[str, NDArray[np.floating]]:
        if not self.points:
            empty = np.array([], dtype=np.float64)
            return {k: empty for k in ("x", "r", "theta", "nu", "M", "mu")}
        return {
            "x": np.array([p.x for p in self.points], dtype=np.float64),
            "r": np.array([p.r for p in self.points], dtype=np.float64),
            "theta": np.array([p.theta for p in self.points], dtype=np.float64),
            "nu": np.array([p.nu for p in self.points], dtype=np.float64),
            "M": np.array([p.M for p in self.points], dtype=np.float64),
            "mu": np.array([p.mu for p in self.points], dtype=np.float64),
        }
