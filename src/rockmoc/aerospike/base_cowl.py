"""Aerospike truncation base pressure + primary cowl (internal nozzle) helpers.

Base pressure on a truncated spike is a major real-world Isp correction.
Models here are engineering correlations, not full CFD.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from rockmoc.aerospike.lee_spike import SpikeResult, design_spike, design_spike_truncated
from rockmoc.exceptions import DesignError
from rockmoc.geometry.contour import Contour
from rockmoc.thermo.perfect_gas import PerfectGas


@dataclass(slots=True)
class BasePressureModel:
    """Simple truncated-spike base pressure estimate."""

    p_base: float
    p_base_over_pc: float
    model: str
    meta: dict

    def summary(self) -> dict:
        return {
            "p_base_Pa": self.p_base,
            "p_base_over_pc": self.p_base_over_pc,
            "model": self.model,
            **self.meta,
        }


def estimate_base_pressure(
    *,
    pc: float,
    pa: float,
    pe: float,
    truncation_fraction: float,
    gamma: float = 1.25,
) -> BasePressureModel:
    """Estimate base pressure on truncated aerospike tip.

    Correlation (engineering):
        p_b / p_a ≈ k1 + k2 * (1 - L/L_full) * (pe/pa)^n
    clamped between a fraction of pe and a multiple of pa.

    For full spike (truncation_fraction=1), base area → 0, p_b unused.
    """
    if pc <= 0:
        raise DesignError("pc must be positive")
    f = float(np.clip(truncation_fraction, 0.05, 1.0))
    if f >= 0.999:
        return BasePressureModel(
            p_base=pa,
            p_base_over_pc=pa / pc,
            model="full_spike_no_base",
            meta={"truncation_fraction": f},
        )

    # Open literature: base pressure often O(pa) with recovery toward pe
    # for highly truncated plugs.
    trunc = 1.0 - f
    # Simple blend: more truncation → p_b closer to pe (recirculation)
    # less truncation → p_b closer to pa
    pe_use = max(pe, 1e-6)
    pa_use = max(pa, 0.0)
    p_b = pa_use * (1.0 - 0.6 * trunc) + pe_use * (0.35 * trunc)
    # Ambient dominance at high altitude
    if pa_use < 0.01 * pe_use:
        p_b = 0.3 * pe_use * trunc + pa_use
    p_b = max(p_b, 0.0)

    return BasePressureModel(
        p_base=float(p_b),
        p_base_over_pc=float(p_b / pc),
        model="blend_pa_pe_truncation",
        meta={
            "truncation_fraction": f,
            "pa_Pa": pa,
            "pe_Pa": pe,
            "gamma": gamma,
        },
    )


def truncated_spike_with_base(
    full: SpikeResult,
    *,
    length_fraction: float = 0.2,
) -> tuple[SpikeResult, BasePressureModel, dict]:
    """Truncate spike and correct thrust/Isp with base pressure force.

    ΔF_base ≈ (p_b - p_a) * A_base
    A_base ≈ π R_tip²  (projected base of truncated spike).
    """
    trunc = design_spike_truncated(full, length_fraction=length_fraction)
    R_tip = float(trunc.R[-1])
    A_base = math.pi * max(R_tip, 0.0) ** 2
    bp = estimate_base_pressure(
        pc=full.Pc,
        pa=full.Pa,
        pe=full.Pe,
        truncation_fraction=length_fraction,
        gamma=full.gamma,
    )
    dF = (bp.p_base - full.Pa) * A_base
    F_corr = trunc.F + dF
    Isp_corr = F_corr / (trunc.m_dot * 9.80665) if trunc.m_dot > 0 else trunc.Isp
    Cf_corr = F_corr / (full.Pc * full.At) if full.At > 0 else trunc.Cf

    perf = {
        "F_uncorrected_N": trunc.F,
        "F_base_N": dF,
        "F_corrected_N": F_corr,
        "Isp_uncorrected_s": trunc.Isp,
        "Isp_corrected_s": Isp_corr,
        "Cf_corrected": Cf_corr,
        "A_base_m2": A_base,
        "R_tip_m": R_tip,
        "base": bp.summary(),
    }
    # Update truncated result fields
    trunc.F = F_corr
    trunc.Isp = Isp_corr
    trunc.Cf = Cf_corr
    trunc.meta = {**trunc.meta, "base_pressure_correction": perf}
    return trunc, bp, perf


@dataclass(slots=True)
class CowlDesign:
    """Simple primary (internal) cowl nozzle feeding the aerospike lip."""

    contour: Contour
    throat_radius: float  # annular gap scale proxy
    exit_radius: float
    length: float
    area_ratio_primary: float
    meta: dict

    def summary(self) -> dict:
        return {
            "length_m": self.length,
            "area_ratio_primary": self.area_ratio_primary,
            "throat_radius_m": self.throat_radius,
            "exit_radius_m": self.exit_radius,
            **self.meta,
        }


def design_primary_cowl(
    *,
    Re_lip: float,
    ht: float,
    delta_rad: float,
    gamma: float = 1.25,
    primary_eps: float = 2.0,
    n_points: int = 80,
) -> CowlDesign:
    """Design a short conical/primary internal nozzle to the shroud lip.

    Geometry (engineering):
    - Throat annulus height ht at angle delta to the axis (Lee/Vernacchia).
    - Short expansion to primary_eps along the shroud side before the lip.
    - Contour returned in (x, r) with lip at r ≈ Re_lip.
    """
    if Re_lip <= 0 or ht <= 0:
        raise DesignError("Re_lip and ht must be positive")
    if primary_eps < 1.0:
        raise DesignError("primary_eps must be ≥ 1")

    gas = PerfectGas(gamma=gamma)
    # Throat at lip-inner side: r_t_outer ≈ Re_lip, gap ht along shroud normal
    # Approximate centerline of annulus
    r_throat_mid = Re_lip - 0.5 * ht * math.sin(delta_rad)
    # Primary exit at lip
    r_exit = Re_lip
    # Length of primary nozzle along axis
    # Area ratio for annular: rough (r_e^2 - r_i^2)/(r_t^2 - r_i_t^2)
    # Use simple axial length from 15° wall
    L = 1.5 * ht / max(math.tan(abs(delta_rad) + 0.15), 0.1)
    L = max(L, 0.5 * ht)

    x = np.linspace(-L, 0.0, n_points)  # lip at x=0
    # Outer shroud wall: straight from throat station to lip
    r_outer = np.linspace(r_throat_mid + 0.5 * ht, r_exit, n_points)
    # Inner wall (spike side at primary): gap grows with primary_eps
    gap = ht * np.sqrt(np.linspace(1.0, primary_eps, n_points))
    r_inner = r_outer - gap * math.sin(max(delta_rad, 0.2))

    # Export outer cowl contour as the primary wall
    th = np.arctan2(np.gradient(r_outer), np.gradient(x))
    M = np.linspace(1.05, float(gas.mach_from_area_ratio(primary_eps, branch="supersonic")), n_points)
    contour = Contour(x=x, r=r_outer, theta=th, M=M)

    return CowlDesign(
        contour=contour,
        throat_radius=ht,  # gap as characteristic scale
        exit_radius=r_exit,
        length=L,
        area_ratio_primary=primary_eps,
        meta={
            "delta_deg": math.degrees(delta_rad),
            "Re_lip_m": Re_lip,
            "ht_m": ht,
            "r_inner_throat": float(r_inner[0]),
            "r_inner_exit": float(r_inner[-1]),
            "model": "short_primary_cowl_engineering",
        },
    )


def design_spike_with_cowl(
    *,
    Pc: float,
    Tc: float,
    gamma: float,
    molar_mass: float,
    Re: float,
    expansion_ratio: float,
    Pa: float = 101325.0,
    n_points: int = 10000,
    primary_eps: float = 2.0,
    truncate: float | None = None,
) -> dict:
    """Full annular stack: primary cowl + spike (+ optional truncation/base)."""
    spike = design_spike(
        Pc=Pc,
        Tc=Tc,
        gamma=gamma,
        molar_mass=molar_mass,
        Re=Re,
        expansion_ratio=expansion_ratio,
        Pa=Pa,
        n_points=n_points,
    )
    cowl = design_primary_cowl(
        Re_lip=Re,
        ht=spike.ht,
        delta_rad=spike.delta_rad,
        gamma=gamma,
        primary_eps=primary_eps,
    )
    out: dict = {
        "spike": spike,
        "cowl": cowl,
        "spike_summary": spike.summary(),
        "cowl_summary": cowl.summary(),
    }
    if truncate is not None:
        trunc, bp, perf = truncated_spike_with_base(spike, length_fraction=truncate)
        out["truncated_spike"] = trunc
        out["base_pressure"] = bp
        out["base_performance"] = perf
    return out
