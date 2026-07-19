"""Design specification and result types."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

from rockmoc.exceptions import DesignError
from rockmoc.geometry.contour import Contour
from rockmoc.moc.types import CharNet
from rockmoc.performance.metrics import PerformanceMetrics
from rockmoc.thermo.perfect_gas import PerfectGas


@dataclass(slots=True)
class DesignSpec:
    """Inputs for a nozzle contour design.

    Provide **exactly one** of ``exit_mach``, ``area_ratio``, or ``pressure_ratio``
    (``pe/pc``) to set the design expansion.

    Parameters
    ----------
    gas:
        Perfect-gas model.
    throat_radius:
        Geometric throat radius Rt [m].
    n_characteristics:
        Characteristic mesh density (fan rays). Larger → smoother, slower.
    throat_mach:
        Slightly supersonic start Mach on the initial data line (> 1).
    throat_rc_over_rt:
        Throat wall radius of curvature / Rt (ideal / TIC circular arc).
    sharp_corner:
        If True (default for MLN), use a sharp throat corner.
    delta:
        Geometry flag: 1 = axisymmetric (default), 0 = planar (internal/testing).
    """

    gas: PerfectGas = field(default_factory=PerfectGas)  # or VariableGammaGas / MOCGas
    throat_radius: float = 0.05
    exit_mach: float | None = None
    area_ratio: float | None = None
    pressure_ratio: float | None = None  # pe/pc
    n_characteristics: int = 160  # fine mesh default for accuracy
    throat_mach: float = 1.05
    throat_rc_over_rt: float = 1.5
    sharp_corner: bool = True
    delta: float = 1.0
    # TIC only
    truncate_fraction: float | None = None

    def resolved_exit_mach(self) -> float:
        """Resolve design exit Mach from the active targeting mode."""
        n_set = sum(
            x is not None
            for x in (self.exit_mach, self.area_ratio, self.pressure_ratio)
        )
        if n_set != 1:
            raise DesignError(
                "Provide exactly one of exit_mach, area_ratio, or pressure_ratio"
            )
        if self.exit_mach is not None:
            Me = float(self.exit_mach)
        elif self.area_ratio is not None:
            if self.area_ratio < 1.0:
                raise DesignError("area_ratio must be ≥ 1")
            Me = float(
                self.gas.mach_from_area_ratio(self.area_ratio, branch="supersonic")
            )
        else:
            assert self.pressure_ratio is not None
            pe_pc = float(self.pressure_ratio)
            if not (0.0 < pe_pc < 1.0):
                raise DesignError("pressure_ratio pe/pc must be in (0, 1)")
            # pe/pc = (1 + 0.5(γ-1)M^2)^(-γ/(γ-1))
            g = self.gas.gamma
            T_ratio = pe_pc ** (-(g - 1.0) / g)
            Me = math_sqrt_mach_from_temp_ratio(T_ratio, g)

        if Me <= 1.0:
            raise DesignError(f"Design exit Mach must be > 1; got {Me}")
        return Me

    def validate(self) -> None:
        if self.throat_radius <= 0.0:
            raise DesignError("throat_radius must be positive")
        if self.n_characteristics < 4:
            raise DesignError("n_characteristics must be ≥ 4")
        if self.throat_mach <= 1.0:
            raise DesignError("throat_mach must be > 1")
        if self.throat_rc_over_rt <= 0.0:
            raise DesignError("throat_rc_over_rt must be positive")
        if self.delta not in (0.0, 1.0):
            raise DesignError("delta must be 0 (planar) or 1 (axisymmetric)")
        if self.truncate_fraction is not None and not (
            0.05 < self.truncate_fraction <= 1.0
        ):
            raise DesignError("truncate_fraction must be in (0.05, 1]")
        _ = self.resolved_exit_mach()


def math_sqrt_mach_from_temp_ratio(T0_over_T: float, gamma: float) -> float:
    """M from T0/T = 1 + 0.5(γ-1)M²."""
    if T0_over_T < 1.0:
        raise DesignError("Invalid temperature ratio for Mach solve")
    return float(np.sqrt(2.0 * (T0_over_T - 1.0) / (gamma - 1.0)))


@dataclass(slots=True)
class DesignResult:
    """Output of a nozzle design."""

    kind: Literal["mln", "ideal", "tic", "rao"]
    spec: DesignSpec
    exit_mach_design: float
    contour: Contour
    metrics: PerformanceMetrics
    net: CharNet | None = None
    meta: dict[str, Any] = field(default_factory=dict)

    def export_csv(self, path: str) -> None:
        self.contour.export_csv(path)

    def export_json(self, path: str) -> None:
        self.contour.export_json(path, extra=self.summary())

    def summary(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "exit_mach_design": self.exit_mach_design,
            "gamma": self.spec.gas.gamma,
            "throat_radius_m": self.spec.throat_radius,
            "n_characteristics": self.spec.n_characteristics,
            "area_ratio": self.metrics.area_ratio,
            "length_m": self.metrics.length,
            "length_over_Rt": self.metrics.length / self.spec.throat_radius,
            "exit_radius_m": self.metrics.exit_radius,
            "Cf_vac_1d": self.metrics.Cf_vac_1d,
            "max_wall_angle_deg": float(np.degrees(self.metrics.max_wall_angle)),
            **self.meta,
        }

    def plot(self, **kwargs: Any) -> Any:
        """Show dark-themed nozzle + MOC dashboard after solve.

        Parameters (kwargs)
        -------------------
        show : bool
            Display figures (default True).
        save_dir : str | None
            If set, save PNGs into this directory (optional).
        prefix : str
            Filename prefix when saving.
        """
        from rockmoc.plot.dashboard import display_design_result

        show = kwargs.pop("show", True)
        save_dir = kwargs.pop("save_dir", None)
        prefix = kwargs.pop("prefix", self.kind)
        return display_design_result(
            self, show=show, save_dir=save_dir, prefix=prefix, **kwargs
        )
