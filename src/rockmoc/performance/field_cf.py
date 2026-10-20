"""Field-integrated thrust coefficient from MOC net / exit plane.

Integrates pressure and axial momentum on a discrete exit ray (axis → wall)
when a characteristic net is available; falls back to contour-edge 1-D +
divergence efficiency when not.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from rockmoc.geometry.contour import Contour
from rockmoc.moc.types import CharNet
from rockmoc.performance.integrated import IntegratedPerformance, divergence_efficiency
from rockmoc.thermo.perfect_gas import PerfectGas

G0 = 9.80665


@dataclass(slots=True)
class FieldCfResult:
    """Exit-plane integrated performance."""

    Cf_vac: float
    Cf_amb: float
    Isp_vac: float
    Isp_amb: float
    pe_pc_mean: float
    Me_mean: float
    Me_axis: float | None
    Me_wall: float | None
    theta_exit_max: float
    mass_flow: float
    thrust_N: float
    method: str
    meta: dict

    def summary(self) -> dict:
        return {
            "method": self.method,
            "Cf_vac": self.Cf_vac,
            "Cf_amb": self.Cf_amb,
            "Isp_vac_s": self.Isp_vac,
            "Isp_amb_s": self.Isp_amb,
            "pe_pc_mean": self.pe_pc_mean,
            "Me_mean": self.Me_mean,
            "Me_axis": self.Me_axis,
            "Me_wall": self.Me_wall,
            "theta_exit_max_deg": math.degrees(self.theta_exit_max),
            "m_dot_kg_s": self.mass_flow,
            "thrust_N": self.thrust_N,
            **self.meta,
        }


def _exit_ray_from_net(net: CharNet, *, n: int = 64) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """(r, M, θ) on the exit plane x = x_lip, interpolated from the characteristic net.

    Linear interpolation on the Delaunay triangulation of all net points
    (axis, interior, transition, wall); the lip is the net's last wall node.
    Returns None for nets that do not reach the exit (e.g. kernel-only).
    """
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

    if not net.wall or len(net.points) < 10:
        return None
    lip = net.wall[-1]
    pts = list({id(p): p for p in [*net.points, *net.wall, *net.axis]}.values())
    xy = np.array([(p.x, p.r) for p in pts])
    if xy[:, 0].max() < lip.x * (1.0 - 1e-6) or not np.any(np.abs(xy[:, 0] - lip.x) < 0.05 * lip.r):
        return None
    vals = np.array([(p.M, p.theta) for p in pts])
    r = np.linspace(0.0, lip.r, n)
    q = np.column_stack([np.full(n, lip.x), r])
    out = LinearNDInterpolator(xy, vals)(q)
    bad = np.isnan(out[:, 0])
    if bad.any():
        out[bad] = NearestNDInterpolator(xy, vals)(q[bad])
    out[-1] = (lip.M, lip.theta)
    out[0, 1] = 0.0
    return r, out[:, 0], out[:, 1]


def _exit_station_from_contour(
    contour: Contour,
    *,
    exit_mach_design: float,
    n: int = 48,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Axisymmetric exit ray: r ∈ [0, R_e], linear θ(r), isentropic M(r).

    M(r) is taken uniform at design Me when contour.M is missing; otherwise
    a mild core-to-wall blend from Me_design to wall M preserves edge state.
    """
    Re = float(contour.r[-1])
    r = np.linspace(0.0, Re, n)
    th_w = float(contour.theta[-1]) if contour.theta is not None else 0.0
    th = th_w * (r / max(Re, 1e-12))
    Me = float(exit_mach_design)
    if contour.M is not None and len(contour.M) > 0:
        Mw = float(contour.M[-1])
        # Slightly higher Mach on axis is typical for contoured bells
        M_axis = max(Me, Mw)
        M = M_axis + (Mw - M_axis) * (r / max(Re, 1e-12)) ** 2
    else:
        M = np.full_like(r, Me)
    return r, M, th, "contour_exit_station"


def field_integrated_cf(
    contour: Contour,
    gas: PerfectGas,
    *,
    throat_radius: float,
    exit_mach_design: float,
    pc: float = 5e6,
    Tc: float = 3000.0,
    pa: float = 0.0,
    net: CharNet | None = None,
) -> FieldCfResult:
    """Integrate thrust on the nozzle exit plane (axisymmetric).

    Exit-plane integrals (absolute pressure form)::

        ṁ = ∫ ρ u  2π r dr
        F  = ∫ (ρ u² + p) 2π r dr  − p_a A_e
        u  = V cos θ

    States are isentropic from chamber (p_c, T_c) using local M. Divergence is
    carried by θ(r) inside the integral — no extra η_div factor is applied.
    """
    Rt = throat_radius
    g = gas.gamma
    At = math.pi * Rt**2

    ray = _exit_ray_from_net(net) if net is not None else None
    if net is not None and ray is None:
        from rockmoc.exceptions import warn_fallback

        warn_fallback("MOC net does not reach the exit plane; synthetic exit station used for Cf")
    if ray is not None:
        method = "moc_exit_plane"
        r, M, th = ray
        # The net may be at the un-rescaled (mass-consistent) size; scale the
        # exit plane to this contour's exit radius so At / Ae stay consistent.
        r = r * (float(contour.r[-1]) / float(r[-1]))
    else:
        r, M, th, method = _exit_station_from_contour(
            contour, exit_mach_design=exit_mach_design
        )

    M = np.maximum(M, 1.01)
    T = Tc / (1.0 + 0.5 * (g - 1.0) * M**2)
    p = pc / (1.0 + 0.5 * (g - 1.0) * M**2) ** (g / (g - 1.0))
    rho = p / (gas.R * T)
    a = np.sqrt(g * gas.R * T)
    V = M * a
    u = V * np.cos(th)

    # Trapezoidal axisymmetric integrals
    integrand_m = 2.0 * math.pi * r * rho * u
    integrand_F = 2.0 * math.pi * r * (rho * u**2 + p)
    trapz = getattr(np, "trapezoid", None) or np.trapz
    m_dot = float(trapz(integrand_m, r))
    F_abs = float(trapz(integrand_F, r))
    Ae = math.pi * float(r[-1]) ** 2
    F = F_abs - pa * Ae
    F_vac = F_abs

    # Throat mass flow for reference (ideal 1-D)
    cstar = gas.characteristic_velocity(Tc)
    m_dot_1d = pc * At / cstar if cstar > 0 else m_dot
    # Prefer consistent Cf from force / (pc At); mass for Isp from integral
    # renormalized to throat continuity when ray m_dot is off (sparse state)
    mass_error = (m_dot - m_dot_1d) / max(m_dot_1d, 1e-12)
    if method != "moc_exit_plane" and m_dot > 1e-12 and abs(mass_error) > 0.15:
        # Scale density level so m_dot matches 1-D throat (same M,θ shape)
        scale_m = m_dot_1d / m_dot
        m_dot = m_dot_1d
        F_abs *= scale_m
        # pressure part of F should not all scale with density; re-integrate
        # split: F = F_mom + F_p − pa Ae
        F_mom = float(trapz(2.0 * math.pi * r * rho * u**2, r)) * scale_m
        F_p = float(trapz(2.0 * math.pi * r * p, r))
        F_abs = F_mom + F_p
        F_vac = F_abs
        F = F_abs - pa * Ae

    if method == "moc_exit_plane" and m_dot > 1e-12:
        # C_F = F / (ṁ c*): uses the field's own mass flow, so the result is
        # independent of the net's absolute scale (discharge coefficient, etc.)
        Cf_vac = F_vac / (m_dot * cstar)
        Cf_amb = F / (m_dot * cstar)
    else:
        Cf_vac = F_vac / (pc * At) if pc * At > 0 else float("nan")
        Cf_amb = F / (pc * At) if pc * At > 0 else float("nan")

    if m_dot > 1e-12:
        Isp_vac = F_vac / (m_dot * G0)
        Isp_amb = F / (m_dot * G0)
    else:
        Isp_vac = Cf_vac * cstar / G0
        Isp_amb = Cf_amb * cstar / G0

    pe_mean = float(np.mean(p))
    return FieldCfResult(
        Cf_vac=float(Cf_vac),
        Cf_amb=float(Cf_amb),
        Isp_vac=float(Isp_vac),
        Isp_amb=float(Isp_amb),
        pe_pc_mean=pe_mean / pc if pc > 0 else 0.0,
        Me_mean=float(np.mean(M)),
        Me_axis=float(M[0]) if len(M) else None,
        Me_wall=float(M[-1]) if len(M) else None,
        theta_exit_max=float(np.max(np.abs(th))) if len(th) else 0.0,
        mass_flow=m_dot,
        thrust_N=float(F),
        method=method,
        meta={
            "pc_Pa": pc,
            "Tc_K": Tc,
            "pa_Pa": pa,
            "At_m2": At,
            "Ae_m2": Ae,
            "cstar_m_s": float(cstar),
            "m_dot_1d_kg_s": float(m_dot_1d),
            "n_ray_points": int(len(r)),
            "mass_error_vs_1d": float(mass_error),
            "exit_mach_design": exit_mach_design,
            "eta_div_contour": divergence_efficiency(contour),
        },
    )


def field_to_integrated(field: FieldCfResult) -> IntegratedPerformance:
    """Adapt FieldCfResult to IntegratedPerformance shape."""
    return IntegratedPerformance(
        Cf_vac=field.Cf_vac,
        Cf_amb=field.Cf_amb,
        Isp_vac=field.Isp_vac,
        Isp_amb=field.Isp_amb,
        pe_pc=field.pe_pc_mean,
        Me=field.Me_mean,
        area_ratio=field.meta.get("Ae_m2", 0.0)
        / max(field.meta.get("At_m2", 1.0), 1e-30)
        if field.meta.get("At_m2")
        else 0.0,
        divergence_efficiency=1.0,
        kinetic_efficiency=1.0,
        nozzle_efficiency=1.0,
        meta={"field": field.summary()},
    )
