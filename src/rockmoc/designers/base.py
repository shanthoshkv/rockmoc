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

    @classmethod
    def student(
        cls,
        *,
        Me: float = 3.0,
        gamma: float = 1.25,
        Rt_mm: float = 50.0,
        R: float = 350.0,
        n: int = 40,
        gas_name: str | None = None,
        **kwargs: Any,
    ) -> DesignSpec:
        """Homework-friendly constructor (mm throat, modest mesh).

        Parameters
        ----------
        Me:
            Design exit Mach.
        gamma:
            Specific-heat ratio (ignored if ``gas_name`` is set).
        Rt_mm:
            Throat radius in **millimeters** (converted to meters).
        R:
            Specific gas constant [J/(kg·K)].
        n:
            Characteristic count (default 40 — interactive speed).
        gas_name:
            Optional preset from :func:`rockmoc.student.example_gas`
            (e.g. ``\"lox_rp1\"``, ``\"air\"``).
        """
        if gas_name is not None:
            from rockmoc.student.gases import example_gas

            gas = example_gas(gas_name)
        else:
            gas = PerfectGas(gamma=gamma, R=R)
        return cls(
            gas=gas,
            throat_radius=float(Rt_mm) * 1e-3,
            exit_mach=float(Me),
            n_characteristics=int(n),
            **kwargs,
        )

    def resolved_exit_mach(self) -> float:
        """Resolve design exit Mach from the active targeting mode."""
        n_set = sum(
            x is not None
            for x in (self.exit_mach, self.area_ratio, self.pressure_ratio)
        )
        if n_set != 1:
            set_names = [
                n
                for n, v in (
                    ("exit_mach", self.exit_mach),
                    ("area_ratio", self.area_ratio),
                    ("pressure_ratio", self.pressure_ratio),
                )
                if v is not None
            ]
            raise DesignError(
                "Provide exactly one of exit_mach, area_ratio, or pressure_ratio "
                f"(you set {set_names or 'none'}).",
                tip="Leave the unused targets as None.",
                example="DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0)",
            )
        if self.exit_mach is not None:
            Me = float(self.exit_mach)
        elif self.area_ratio is not None:
            if self.area_ratio < 1.0:
                raise DesignError(
                    f"area_ratio must be ≥ 1; got {self.area_ratio}",
                    tip="ε = Ae/A* is at least 1 at the throat (M=1).",
                    example="DesignSpec(..., area_ratio=25.0)",
                )
            Me = float(
                self.gas.mach_from_area_ratio(self.area_ratio, branch="supersonic")
            )
        else:
            assert self.pressure_ratio is not None
            pe_pc = float(self.pressure_ratio)
            if not (0.0 < pe_pc < 1.0):
                raise DesignError(
                    f"pressure_ratio pe/pc must be in (0, 1); got {pe_pc}",
                    tip="Use exit static pressure over chamber pressure (e.g. 0.01).",
                    example="DesignSpec(..., pressure_ratio=0.01)",
                )
            # pe/pc = (1 + 0.5(γ-1)M^2)^(-γ/(γ-1))
            g = self.gas.gamma
            T_ratio = pe_pc ** (-(g - 1.0) / g)
            Me = math_sqrt_mach_from_temp_ratio(T_ratio, g)

        if Me <= 1.0:
            raise DesignError(
                f"Design exit Mach must be > 1; got {Me}",
                tip="Supersonic nozzles need Me > 1 (try 2–4 for homework).",
                example="DesignSpec(..., exit_mach=3.0)",
            )
        return Me

    def validate(self) -> None:
        if self.throat_radius <= 0.0:
            raise DesignError(
                f"throat_radius must be positive; got {self.throat_radius}",
                tip="Use meters (e.g. 0.05 for 50 mm), or DesignSpec.student(Rt_mm=50).",
                example="DesignSpec(..., throat_radius=0.05)",
            )
        if self.n_characteristics < 4:
            raise DesignError(
                f"n_characteristics must be ≥ 4; got {self.n_characteristics}",
                tip="Use 24–40 for quick demos, 80–160 for finer walls.",
                example="DesignSpec(..., n_characteristics=40)",
            )
        if self.throat_mach <= 1.0:
            raise DesignError(
                f"throat_mach must be > 1; got {self.throat_mach}",
                tip="Start slightly supersonic (default 1.05).",
            )
        if self.throat_rc_over_rt <= 0.0:
            raise DesignError(
                "throat_rc_over_rt must be positive",
                tip="Typical values are 0.5–2.0 (default 1.5).",
            )
        if self.delta not in (0.0, 1.0):
            raise DesignError(
                "delta must be 0 (planar) or 1 (axisymmetric)",
                tip="Leave delta=1.0 for normal rocket bells.",
            )
        if self.truncate_fraction is not None and not (
            0.05 < self.truncate_fraction <= 1.0
        ):
            raise DesignError(
                f"truncate_fraction must be in (0.05, 1]; got {self.truncate_fraction}",
                tip="For TIC try 0.8 (80% of ideal length).",
                example="DesignSpec(..., truncate_fraction=0.8)",
            )
        _ = self.resolved_exit_mach()


def math_sqrt_mach_from_temp_ratio(T0_over_T: float, gamma: float) -> float:
    """M from T0/T = 1 + 0.5(γ-1)M²."""
    if T0_over_T < 1.0:
        raise DesignError("Invalid temperature ratio for Mach solve")
    return float(np.sqrt(2.0 * (T0_over_T - 1.0) / (gamma - 1.0)))


@dataclass(slots=True)
class DesignResult:
    """Output of a nozzle design."""

    kind: Literal["mln", "ideal", "tic", "rao", "conical"]
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
        from rockmoc.student.limits import LIMITS_SHORT

        th_e = (
            float(np.degrees(abs(self.contour.theta[-1])))
            if self.contour.theta is not None
            else 0.0
        )
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
            "exit_wall_angle_deg": th_e,
            "limitations": LIMITS_SHORT,
            **self.meta,
        }

    def report(
        self,
        *,
        length_unit: str = "mm",
        angle_unit: str = "deg",
        include_limits: bool = True,
        include_explainer: bool = True,
    ) -> str:
        """Unit-aware lab-style text report (print this in notebooks)."""
        from rockmoc.student.report import format_design_report

        return format_design_report(
            self,
            length_unit=length_unit,  # type: ignore[arg-type]
            angle_unit=angle_unit,  # type: ignore[arg-type]
            include_limits=include_limits,
            include_explainer=include_explainer,
        )

    def lab_report(
        self,
        path: str,
        *,
        student_name: str = "",
        course: str = "",
        length_unit: str = "mm",
        angle_unit: str = "deg",
        extra_notes: str = "",
    ) -> str:
        """Write a Markdown lab-report skeleton with checklist; return path."""
        from rockmoc.student.report import generate_lab_report

        return generate_lab_report(
            self,
            path,
            student_name=student_name,
            course=course,
            length_unit=length_unit,  # type: ignore[arg-type]
            angle_unit=angle_unit,  # type: ignore[arg-type]
            extra_notes=extra_notes,
        )

    def with_boundary_layer(
        self, *, pc: float = 5.0e6, Tc: float = 3000.0, turbulent: bool = True, **bl_kwargs: Any
    ) -> DesignResult:
        """Copy of this design with the wall offset outward by the displacement thickness δ*.

        Uses :func:`rockmoc.analysis.correct_wall_bl` (momentum-integral BL
        with empirical fallback); extra keywords (``method``,
        ``wall_temperature``, …) are forwarded. δ*, friction ΔCf and the
        method actually used are added to ``meta["boundary_layer"]``.
        """
        from rockmoc.analysis.boundary_layer import correct_wall_bl
        from rockmoc.performance.metrics import compute_metrics

        bl = correct_wall_bl(
            self.contour, self.spec.gas, throat_radius=self.spec.throat_radius,
            pc=pc, Tc=Tc, turbulent=turbulent, axisymmetric=self.spec.delta == 1.0,
            **bl_kwargs,
        )
        metrics = compute_metrics(
            bl.corrected, self.spec.gas, throat_radius=self.spec.throat_radius,
            exit_mach_design=self.exit_mach_design,
        )
        return DesignResult(
            kind=self.kind, spec=self.spec, exit_mach_design=self.exit_mach_design,
            contour=bl.corrected, metrics=metrics, net=self.net,
            meta={**self.meta, "boundary_layer": bl.summary()},
        )

    def explain(self) -> str:
        """Printable mini-explainer for this nozzle kind."""
        from rockmoc.student.explainers import explainer

        return explainer(str(self.kind))

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
