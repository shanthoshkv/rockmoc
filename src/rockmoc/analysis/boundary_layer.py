"""Boundary-layer displacement correction for inviscid MOC walls.

Primary reference (local copy ``refs/NASA-TM-X-2343.pdf``)
----------------------------------------------------------
Goldman & Vanco, NASA TM-X-2343 (1971): design sharp-edged-throat nozzles by
MOC, then **add displacement thickness** to ideal coordinates. Laminar:
Cohen–Reshotko; turbulent: Sasman–Cresci. Curvature neglected.

This module implements a compressible integral estimate (not the full
Cohen–Reshotko / Sasman–Cresci ODEs) that follows the same *geometry*
rule: r_corrected = r_ideal + δ* (along the wall normal / radial).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from rockmoc.exceptions import DesignError
from rockmoc.geometry.contour import Contour
from rockmoc.thermo.perfect_gas import PerfectGas


@dataclass(slots=True)
class BLCorrectionResult:
    """Inviscid + BL-corrected contours and δ* distribution."""

    inviscid: Contour
    corrected: Contour
    delta_star: NDArray[np.floating]  # m
    theta_mom: NDArray[np.floating]  # momentum thickness m
    Re_x: NDArray[np.floating]
    meta: dict

    def summary(self) -> dict:
        return {
            "method": "integral_bl_displacement",
            "delta_star_exit_m": float(self.delta_star[-1]),
            "delta_star_max_m": float(np.max(self.delta_star)),
            "exit_radius_inviscid_m": float(self.inviscid.r[-1]),
            "exit_radius_corrected_m": float(self.corrected.r[-1]),
            "length_m": float(self.corrected.length),
            **self.meta,
        }


def _viscosity_sutherland(T: float, mu0: float = 1.716e-5, T0: float = 273.15, S: float = 110.4) -> float:
    """Air-like Sutherland law [Pa·s]. Acceptable engineering default for hot gas order-of-magnitude."""
    return mu0 * (T / T0) ** 1.5 * (T0 + S) / (T + S)


def correct_wall_bl(
    contour: Contour,
    gas: PerfectGas,
    *,
    throat_radius: float,
    pc: float = 5.0e6,
    Tc: float = 3000.0,
    prandtl: float = 0.72,
    turbulent: bool = True,
    axisymmetric: bool = True,
) -> BLCorrectionResult:
    """Apply displacement-thickness correction to an inviscid wall contour.

    Parameters
    ----------
    contour:
        Inviscid MOC wall (x along axis, r radial).
    gas:
        Perfect-gas model (γ, R).
    pc, Tc:
        Chamber stagnation p [Pa], T [K] for free-stream estimates.
    turbulent:
        If True use 1/7-power compressible estimate; else laminar Blasius-like.
    axisymmetric:
        Grow radius by δ* (axisym wall). For planar, grow half-height by δ*.
    """
    if contour.x.size < 3:
        raise DesignError("contour too short for BL correction")
    if pc <= 0 or Tc <= 0:
        raise DesignError("pc and Tc must be positive")

    x = np.asarray(contour.x, dtype=np.float64)
    r = np.asarray(contour.r, dtype=np.float64)
    n = len(x)

    # Local edge Mach from contour.M if present, else from area ratio
    if contour.M is not None and len(contour.M) == n:
        M_e = np.maximum(np.asarray(contour.M, dtype=np.float64), 1.01)
    else:
        M_e = np.ones(n)
        for i in range(n):
            ar = max((r[i] / throat_radius) ** 2 if axisymmetric else (r[i] / throat_radius), 1.0)
            try:
                M_e[i] = float(gas.mach_from_area_ratio(ar, branch="supersonic"))
            except Exception:
                M_e[i] = 1.01 + 2.0 * i / max(n - 1, 1)

    # Isentropic edge T, p, ρ, U
    g = gas.gamma
    T_e = Tc / (1.0 + 0.5 * (g - 1.0) * M_e**2)
    p_e = pc / (1.0 + 0.5 * (g - 1.0) * M_e**2) ** (g / (g - 1.0))
    rho_e = p_e / (gas.R * T_e)
    a_e = np.sqrt(g * gas.R * T_e)
    U_e = M_e * a_e

    # Arc length from throat
    ds = np.hypot(np.diff(x), np.diff(r))
    s = np.concatenate([[0.0], np.cumsum(ds)])
    s = np.maximum(s, 1e-9)

    mu_e = np.array([_viscosity_sutherland(float(T)) for T in T_e])
    Re_x = rho_e * U_e * s / mu_e
    Re_x = np.maximum(Re_x, 1.0)

    # Recovery / wall temperature estimate (adiabatic wall)
    r_rec = prandtl ** (1.0 / 3.0) if turbulent else math.sqrt(prandtl)
    T_aw = T_e * (1.0 + r_rec * 0.5 * (g - 1.0) * M_e**2)
    # Reference temperature (Eckert)
    T_ref = 0.5 * (T_e + T_aw) + 0.22 * (T_aw - T_e)
    mu_ref = np.array([_viscosity_sutherland(float(T)) for T in T_ref])
    rho_ref = p_e / (gas.R * T_ref)
    Re_ref = np.maximum(rho_ref * U_e * s / mu_ref, 1.0)

    # Momentum thickness θ
    if turbulent:
        # Compressible 1/7-power: θ/x ≈ 0.036 Re_x^{-1/5} * (T_ref/T_e) factors
        cf = 0.0592 * Re_ref ** (-0.2)
        # von Karman: dθ/dx + … ≈ cf/2; integrate approximately θ ≈ 0.036 x Re^{-1/5}
        theta = 0.036 * s * Re_ref ** (-0.2)
        # Shape factor H = δ*/θ ~ 1.3–1.4 incompressible; compressibility increases H
        H = 1.4 * (1.0 + 0.1 * M_e**2)
    else:
        theta = 0.664 * s * Re_ref ** (-0.5)
        H = 2.6 * (1.0 + 0.08 * M_e**2)

    delta_star = H * theta

    # Grow wall along local outward normal (approx radial for small wall angles)
    if contour.theta is not None:
        alpha = np.asarray(contour.theta, dtype=np.float64)
    else:
        alpha = np.arctan2(np.gradient(r), np.gradient(x))

    # Outward normal in (x,r): for wall tangent (1, tan α), normal ~ (-sin α, cos α)
    # Displacement applied primarily in r for nozzle walls
    r_corr = r + delta_star * np.cos(alpha)
    x_corr = x - delta_star * np.sin(alpha)
    # Ensure monotonic x
    for i in range(1, n):
        if x_corr[i] <= x_corr[i - 1]:
            x_corr[i] = x_corr[i - 1] + 1e-12
    r_corr = np.maximum(r_corr, throat_radius * 0.99)

    corrected = Contour(x=x_corr, r=r_corr, theta=alpha, M=M_e)

    return BLCorrectionResult(
        inviscid=contour,
        corrected=corrected,
        delta_star=delta_star,
        theta_mom=theta,
        Re_x=Re_x,
        meta={
            "turbulent": turbulent,
            "pc_Pa": pc,
            "Tc_K": Tc,
            "prandtl": prandtl,
            "axisymmetric": axisymmetric,
            "reference": "integral BL / TM-X-2343 spirit",
        },
    )
