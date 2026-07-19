"""Performance metrics from designed contours."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rockmoc.geometry.contour import Contour
from rockmoc.thermo.perfect_gas import PerfectGas


@dataclass(slots=True)
class PerformanceMetrics:
    """Summary performance and geometry metrics for a designed nozzle."""

    area_ratio: float
    length: float
    exit_radius: float
    throat_radius: float
    max_wall_angle: float  # rad
    exit_mach_1d: float
    pe_pc_1d: float
    Cf_vac_1d: float
    Cf_sl_1d: float | None = None
    surface_area: float | None = None  # revolved wall surface [m²]

    def as_dict(self) -> dict[str, float | None]:
        return {
            "area_ratio": self.area_ratio,
            "length_m": self.length,
            "exit_radius_m": self.exit_radius,
            "throat_radius_m": self.throat_radius,
            "max_wall_angle_rad": self.max_wall_angle,
            "exit_mach_1d": self.exit_mach_1d,
            "pe_pc_1d": self.pe_pc_1d,
            "Cf_vac_1d": self.Cf_vac_1d,
            "Cf_sl_1d": self.Cf_sl_1d,
            "surface_area_m2": self.surface_area,
        }


def compute_metrics(
    contour: Contour,
    gas: PerfectGas,
    *,
    throat_radius: float,
    exit_mach_design: float,
    pa_pc: float = 0.0,
) -> PerformanceMetrics:
    """Compute geometry + ideal 1D isentropic performance metrics."""
    Rt = throat_radius
    Re = contour.exit_radius
    eps = (Re / Rt) ** 2
    # Prefer design Mach for thermo; also report area-implied Mach
    Me = exit_mach_design
    pe_pc = float(gas.pe_pc_from_mach(Me))
    Cf_vac = gas.thrust_coefficient(pe_pc, 0.0, exit_mach=Me, area_ratio=eps)
    Cf_sl = None
    if pa_pc > 0.0:
        Cf_sl = gas.thrust_coefficient(pe_pc, pa_pc, exit_mach=Me, area_ratio=eps)

    # Revolved lateral surface area: ∫ 2π r ds
    dx = np.diff(contour.x)
    dr = np.diff(contour.r)
    ds = np.hypot(dx, dr)
    r_mid = 0.5 * (contour.r[:-1] + contour.r[1:])
    surface = float(np.sum(2.0 * np.pi * r_mid * ds))

    max_angle = float(np.max(np.abs(contour.theta))) if contour.theta is not None else 0.0

    return PerformanceMetrics(
        area_ratio=float(eps),
        length=contour.length,
        exit_radius=float(Re),
        throat_radius=float(Rt),
        max_wall_angle=max_angle,
        exit_mach_1d=float(Me),
        pe_pc_1d=pe_pc,
        Cf_vac_1d=float(Cf_vac),
        Cf_sl_1d=float(Cf_sl) if Cf_sl is not None else None,
        surface_area=surface,
    )
