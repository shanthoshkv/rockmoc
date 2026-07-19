"""MOC-aware / contour-integrated performance and off-design analysis."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from rockmoc.geometry.contour import Contour
from rockmoc.thermo.perfect_gas import PerfectGas


G0 = 9.80665


@dataclass(slots=True)
class IntegratedPerformance:
    """Thrust coefficient and efficiency estimates from contour + 1-D edge."""

    Cf_vac: float
    Cf_amb: float
    Isp_vac: float | None
    Isp_amb: float | None
    pe_pc: float
    Me: float
    area_ratio: float
    divergence_efficiency: float
    kinetic_efficiency: float
    nozzle_efficiency: float
    meta: dict

    def summary(self) -> dict:
        return {
            "Cf_vac": self.Cf_vac,
            "Cf_amb": self.Cf_amb,
            "Isp_vac_s": self.Isp_vac,
            "Isp_amb_s": self.Isp_amb,
            "pe_pc": self.pe_pc,
            "Me": self.Me,
            "area_ratio": self.area_ratio,
            "divergence_efficiency": self.divergence_efficiency,
            "kinetic_efficiency": self.kinetic_efficiency,
            "nozzle_efficiency": self.nozzle_efficiency,
            **self.meta,
        }


def divergence_efficiency(contour: Contour) -> float:
    """Exit-plane divergence efficiency from exit wall angle.

    Standard conical nozzle correction (also used for contoured bells when
    only θ_e is known)::

        η_div = (1 + cos θ_e) / 2

    For a well-designed MLN / ideal bell, θ_e ≈ 0 ⇒ η_div ≈ 1.
    """
    if contour.theta is None or len(contour.theta) < 1:
        return 1.0
    th_e = abs(float(contour.theta[-1]))
    return float(0.5 * (1.0 + math.cos(th_e)))


def exit_kinetic_efficiency(contour: Contour, *, n_quad: int = 32) -> float:
    """Approximate axial kinetic efficiency at the exit plane.

    Assumes a linear θ(r) from 0 on axis to θ_e at the wall and uniform |V|,
    then η_ke = (∫ u² r dr) / (∫ V² r dr) with u = V cos θ.

    For θ_e = 0 this is exactly 1 (no artificial wall-history penalty).
    """
    if contour.theta is None or len(contour.theta) < 1:
        return 1.0
    th_e = abs(float(contour.theta[-1]))
    if th_e < 1e-12:
        return 1.0
    Re = float(contour.r[-1])
    if Re <= 0.0:
        return 1.0
    # ∫_0^{Re} cos²(θ_e r/Re) r dr  /  ∫ r dr
    r = np.linspace(0.0, Re, n_quad)
    th = th_e * (r / Re)
    trapz = getattr(np, "trapezoid", None) or np.trapz
    num = float(trapz(np.cos(th) ** 2 * r, r))
    den = float(trapz(r, r))
    if den <= 0.0:
        return 1.0
    return float(np.clip(num / den, 0.90, 1.0))


def integrated_performance(
    contour: Contour,
    gas: PerfectGas,
    *,
    throat_radius: float,
    exit_mach: float,
    pc: float = 5e6,
    Tc: float = 3000.0,
    pa: float = 0.0,
) -> IntegratedPerformance:
    """Combine 1-D isentropic Cf with exit-plane divergence / kinetic factors.

    η_n = η_div(θ_e) · η_ke(θ(r) at exit). Wall-angle history is *not* used
    as a loss proxy (that double-counted turning for contoured bells).
    """
    Rt = throat_radius
    Re = float(contour.r[-1])
    eps = (Re / Rt) ** 2
    Me = float(exit_mach)
    pe_pc = float(gas.pe_pc_from_mach(Me))
    pa_pc = pa / pc if pc > 0 else 0.0

    Cf_1d_vac = gas.thrust_coefficient(pe_pc, 0.0, exit_mach=Me, area_ratio=eps)
    Cf_1d_amb = gas.thrust_coefficient(pe_pc, pa_pc, exit_mach=Me, area_ratio=eps)

    eta_div = divergence_efficiency(contour)
    eta_ke = exit_kinetic_efficiency(contour)
    eta_n = eta_div * eta_ke

    Cf_vac = Cf_1d_vac * eta_n
    Cf_amb = Cf_1d_amb * eta_n

    cstar = gas.characteristic_velocity(Tc)
    Isp_vac = Cf_vac * cstar / G0
    Isp_amb = Cf_amb * cstar / G0

    return IntegratedPerformance(
        Cf_vac=float(Cf_vac),
        Cf_amb=float(Cf_amb),
        Isp_vac=float(Isp_vac),
        Isp_amb=float(Isp_amb),
        pe_pc=pe_pc,
        Me=Me,
        area_ratio=float(eps),
        divergence_efficiency=eta_div,
        kinetic_efficiency=eta_ke,
        nozzle_efficiency=eta_n,
        meta={
            "Cf_1d_vac": float(Cf_1d_vac),
            "cstar_m_s": float(cstar),
            "pc_Pa": pc,
            "Tc_K": Tc,
            "pa_Pa": pa,
            "theta_exit_deg": float(
                math.degrees(abs(float(contour.theta[-1]))) if contour.theta is not None else 0.0
            ),
        },
    )


def off_design_altitude_sweep(
    contour: Contour,
    gas: PerfectGas,
    *,
    throat_radius: float,
    exit_mach: float,
    pc: float = 5e6,
    Tc: float = 3000.0,
    altitudes_m: list[float] | None = None,
) -> list[dict]:
    """Cf / Isp vs altitude using simple atmosphere Pa(h)."""
    from rockmoc.analysis.atmosphere import pressure_at_altitude

    if altitudes_m is None:
        altitudes_m = [0.0, 5e3, 10e3, 15e3, 20e3, 30e3, 50e3, 80e3]
    rows = []
    for h in altitudes_m:
        pa = pressure_at_altitude(h)
        perf = integrated_performance(
            contour, gas, throat_radius=throat_radius, exit_mach=exit_mach, pc=pc, Tc=Tc, pa=pa
        )
        rows.append(
            {
                "altitude_m": h,
                "pa_Pa": pa,
                "Cf": perf.Cf_amb,
                "Isp_s": perf.Isp_amb,
                "pe_pc": perf.pe_pc,
                "pa_pc": pa / pc,
            }
        )
    return rows
