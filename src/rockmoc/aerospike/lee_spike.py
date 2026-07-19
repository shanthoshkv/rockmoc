"""Annular aerospike (external-expansion plug) contour — Lee / Vernacchia.

References
----------
1. C.C. Lee, "FORTRAN Programs for Plug Nozzle Design," MSFC / Brown
   Engineering, March 1963.
2. M. Vernacchia, "Spike Contour Algorithm," Pyralis notes + open code
   https://github.com/mvernacc/aerospike-nozzle-design-gui
3. Local copy: ``refs/Spike_Contour_Algorithm.pdf``

Assumptions (as in the reference algorithm)
-------------------------------------------
* Constant γ and molar mass (perfect gas); adiabatic isentropic expansion.
* Expansion fan from the shroud lip; mass continuity on each Mach-line frustum
  fixes the spike radius R_x(M) [Lee eq. (26)].
* Cylindrical (annular) aerospike — not linear.

Default mesh: N = 10_000 points (matches open reference code; dense enough
for fabrication-level resolution on small engines).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from rockmoc.exceptions import DesignError
from rockmoc.geometry.contour import Contour
from rockmoc.thermo.perfect_gas import PerfectGas

G0 = 9.80665  # m/s²
R_UNIV = 8314.462618  # J/(kmol·K)


def mach_from_area_ratio_stodola(eps: float, gamma: float) -> float:
    """Supersonic Mach from expansion ratio (robust root; Stodola inverse)."""
    gas = PerfectGas(gamma=gamma)
    return float(gas.mach_from_area_ratio(eps, branch="supersonic"))


def majdalani_mach_from_er(eps: float, gamma: float) -> float:
    """Explicit series inversion of Stodola's equation (Majdalani & Maickie).

    Used by the Vernacchia reference code; kept for cross-check against
    the general root finder.
    """
    if eps < 1.0:
        raise DesignError("expansion ratio must be ≥ 1")
    if abs(eps - 1.0) < 1e-14:
        return 1.0
    n = 5
    e = 1.0 / float(eps)
    y = float(gamma)
    B = (y + 1.0) / (y - 1.0)
    k = math.sqrt(0.5 * (y - 1.0))
    u = e ** (1.0 / B) / math.sqrt(1.0 + k * k)
    X = np.zeros(n)
    M = np.zeros(n)
    X[0] = (u * k) ** (B / (1.0 - B))
    M[0] = X[0]
    for i in range(1, n):
        lamb = 1.0 / (
            2.0 * M[i - 1] ** (2.0 / B) * (B - 2.0)
            + M[i - 1] ** 2 * B**2 * k**2 * u**2
        )
        disc = (
            M[i - 1] ** (2.0 + 2.0 / B) * k**2 * u**2 * (B**2 - 4.0 * B + 4.0)
            - M[i - 1] ** 2 * B**2 * k**2 * u**4
            + M[i - 1] ** (4.0 / B) * (2.0 * B - 3.0)
            + 2.0 * M[i - 1] ** (2.0 / B) * u**2 * (2.0 - B)
        )
        if disc < 0.0:
            # fall back
            return mach_from_area_ratio_stodola(eps, gamma)
        X[i] = (
            lamb
            * M[i - 1]
            * B
            * (
                M[i - 1] ** (2.0 / B)
                - M[i - 1] ** 2 * B * k**2 * u**2
                + math.sqrt(disc)
            )
        )
        M[i] = M[i - 1] + X[i]
    Me = float(np.real(M[n - 1]))
    if Me < 1.0 or not math.isfinite(Me):
        return mach_from_area_ratio_stodola(eps, gamma)
    return Me


@dataclass(slots=True)
class SpikeResult:
    """Aerospike design output (Lee / Vernacchia)."""

    gamma: float
    molar_mass: float
    Pc: float
    Tc: float
    Pa: float
    Re: float
    expansion_ratio: float
    exit_mach: float
    delta_rad: float
    ht: float
    At: float
    m_dot: float
    F: float
    Isp: float
    Cf: float
    Pe: float
    # Spike contour in physical units (axial X from lip, radius R)
    X: NDArray[np.floating]
    R: NDArray[np.floating]
    M: NDArray[np.floating]
    P: NDArray[np.floating]
    T: NDArray[np.floating]
    Isp_cum: NDArray[np.floating]
    contour: Contour
    meta: dict = field(default_factory=dict)

    def summary(self) -> dict:
        return {
            "method": "Lee_Vernacchia_annular_aerospike",
            "gamma": self.gamma,
            "molar_mass_g_per_mol": self.molar_mass,
            "Pc_Pa": self.Pc,
            "Tc_K": self.Tc,
            "Pa_Pa": self.Pa,
            "Re_m": self.Re,
            "expansion_ratio": self.expansion_ratio,
            "exit_mach": self.exit_mach,
            "delta_deg": math.degrees(self.delta_rad),
            "ht_m": self.ht,
            "At_m2": self.At,
            "m_dot_kg_s": self.m_dot,
            "F_N": self.F,
            "Isp_s": self.Isp,
            "Cf": self.Cf,
            "Pe_Pa": self.Pe,
            "n_points": int(self.X.size),
            "spike_length_m": float(self.X[-1] - self.X[0]) if self.X.size else 0.0,
            "reference": "C.C. Lee 1963; Vernacchia Spike Contour Algorithm",
            **self.meta,
        }


def design_spike(
    *,
    Pc: float,
    Tc: float,
    gamma: float,
    molar_mass: float,
    Re: float,
    expansion_ratio: float,
    Pa: float = 101325.0,
    n_points: int = 20_000,  # fine mesh default (accuracy-first)
    use_majdalani: bool = True,
) -> SpikeResult:
    """Design an annular aerospike spike contour (Lee / Vernacchia algorithm).

    Parameters
    ----------
    Pc, Tc:
        Chamber pressure [Pa] and temperature [K].
    gamma:
        Ratio of specific heats (constant along nozzle in this model).
    molar_mass:
        Average molar mass [g/mol] (same convention as Vernacchia code).
    Re:
        Shroud lip radius [m].
    expansion_ratio:
        ε = A_e / A_t.
    Pa:
        Ambient pressure [Pa] for performance.
    n_points:
        Number of Mach samples along the spike (default 10000).
    use_majdalani:
        If True, use Majdalani series for M_e; else PerfectGas root solve.
    """
    if Pc <= 0 or Tc <= 0 or Re <= 0:
        raise DesignError("Pc, Tc, Re must be positive")
    if expansion_ratio < 1.0:
        raise DesignError("expansion_ratio must be ≥ 1")
    if n_points < 50:
        raise DesignError("n_points must be ≥ 50 for a meaningful contour")
    if molar_mass <= 0:
        raise DesignError("molar_mass must be positive")

    y = float(gamma)
    gas = PerfectGas(gamma=y)
    er = float(expansion_ratio)

    # Exit Mach
    if use_majdalani:
        Me = majdalani_mach_from_er(er, y)
    else:
        Me = mach_from_area_ratio_stodola(er, y)

    # Exit pressure (isentropic)
    Pe = Pc * (1.0 + 0.5 * (y - 1.0) * Me**2) ** (-y / (y - 1.0))

    # Total turning and shroud angle: δ = π/2 − ν_e
    nu_e = float(gas.prandtl_meyer(Me))
    # Shroud angle δ = π/2 − ν_e requires ν_e < π/2 (Lee / Vernacchia).
    # High Me and low γ can violate this ideal full-spike assumption.
    if nu_e >= 0.5 * math.pi - 1e-9:
        raise DesignError(
            f"Exit Prandtl–Meyer angle ν_e={nu_e:.4f} rad ≥ π/2; the Lee/"
            f"Vernacchia ideal spike model requires ν_e < π/2 (reduce "
            f"expansion ratio or increase γ). Me={Me:.3f}, γ={y}."
        )
    delta = 0.5 * math.pi - nu_e
    if delta <= 0.0:
        raise DesignError(
            f"Non-physical shroud angle δ={delta:.4f} rad; check Me/γ "
            f"(ν_e={nu_e:.4f})"
        )

    # Throat gap / Re (Lee / Vernacchia)
    # h_t/R_e = (ε − √[ε(ε − sin δ)]) / (ε sin δ)
    sdel = math.sin(delta)
    if abs(sdel) < 1e-14:
        raise DesignError("sin(δ) ~ 0; cannot form throat gap")
    disc = er * (er - sdel)
    if disc < 0.0:
        raise DesignError("invalid (ε, δ) combination for throat gap formula")
    ht_over_Re = (er - math.sqrt(disc)) / (er * sdel)
    ht = ht_over_Re * Re

    # Specific gas constant [J/(kg·K)]; molar_mass in g/mol
    R = R_UNIV / molar_mass

    # Throat thermo (M=1)
    Tt = Tc / (1.0 + 0.5 * (y - 1.0))
    Pt = Pc * (1.0 + 0.5 * (y - 1.0)) ** (-y / (y - 1.0))
    vt = math.sqrt(y * R * Tt)
    rho_t = Pt / (R * Tt)

    # Throat area (annular): A_t = π h_t (2 R_e − h_t sin δ)
    At = math.pi * ht * (2.0 * Re - ht * sdel)
    m_dot = rho_t * At * vt

    # Mach sweep
    N = int(n_points)
    M = np.linspace(1.0, Me, N)
    mu = np.arcsin(1.0 / M)
    # Prandtl–Meyer (vectorized via gas helper)
    nu = np.array([float(gas.prandtl_meyer(m)) for m in M], dtype=np.float64)
    phi = nu_e - nu + mu

    # Lee eq. (26) / Vernacchia:
    # Rx/Re = sqrt( max(0, 1 − (2/(γ+1)*(1+(γ−1)/2 M²))^((γ+1)/(2(γ−1)))
    #                 * sin(φ) / ε ) )
    expo = (y + 1.0) / (2.0 * (y - 1.0))
    term = (2.0 / (y + 1.0) * (1.0 + 0.5 * (y - 1.0) * M**2)) ** expo
    RxRe2 = 1.0 - term * np.sin(phi) / er
    RxRe2 = np.maximum(RxRe2, 0.0)
    Rx_over_Re = np.sqrt(RxRe2)

    # X/Re = (1 − Rx/Re) / tan(φ)
    tan_phi = np.tan(phi)
    # Avoid division by zero at first point
    tan_phi = np.where(np.abs(tan_phi) < 1e-15, np.sign(tan_phi) * 1e-15, tan_phi)
    X_over_Re = (1.0 - Rx_over_Re) / tan_phi

    R = Rx_over_Re * Re
    X = X_over_Re * Re

    # Thermo along spike
    P = Pc * (1.0 + 0.5 * (y - 1.0) * M**2) ** (-y / (y - 1.0))
    T = Tc / (1.0 + 0.5 * (y - 1.0) * M**2)

    # Cumulative Isp (Vernacchia / Lee finite difference)
    Isp = np.zeros(N, dtype=np.float64)
    Isp[0] = (
        vt
        * math.sin(delta)
        * (
            1.0
            + (1.0 / y)
            * (
                1.0
                - ((y + 1.0) / 2.0) ** (y / (y - 1.0)) * Pa / Pc
            )
        )
        / G0
    )
    pref = vt / y * ((y + 1.0) / 2.0) ** (y / (y - 1.0)) * er / 2.0
    for x in range(1, N):
        Isp[x] = Isp[x - 1] + (
            pref
            * ((P[x - 1] - Pa) / Pc + (P[x] - Pa) / Pc)
            * (Rx_over_Re[x - 1] ** 2 - Rx_over_Re[x] ** 2)
        ) / G0

    F = float(Isp[-1] * m_dot * G0)
    Cf = F / (Pc * At) if At > 0 else float("nan")

    # Contour: axial X as x, radius R as r
    # Ensure finite and sorted
    order = np.argsort(X, kind="mergesort")
    Xo, Ro, Mo = X[order], R[order], M[order]
    # drop non-positive radii at the end (tip)
    contour = Contour(x=Xo, r=np.maximum(Ro, 0.0), theta=None, M=Mo)

    return SpikeResult(
        gamma=y,
        molar_mass=float(molar_mass),
        Pc=float(Pc),
        Tc=float(Tc),
        Pa=float(Pa),
        Re=float(Re),
        expansion_ratio=er,
        exit_mach=float(Me),
        delta_rad=float(delta),
        ht=float(ht),
        At=float(At),
        m_dot=float(m_dot),
        F=float(F),
        Isp=float(Isp[-1]),
        Cf=float(Cf),
        Pe=float(Pe),
        X=Xo,
        R=Ro,
        M=Mo,
        P=P[order],
        T=T[order],
        Isp_cum=Isp[order],
        contour=contour,
        meta={"n_points_requested": N, "use_majdalani": use_majdalani},
    )


def design_spike_truncated(
    result: SpikeResult,
    *,
    length_fraction: float = 0.2,
) -> SpikeResult:
    """Truncate an existing spike design to a fraction of full length.

    Recomputes F/Isp from the cumulative Isp array at the truncate station.
    """
    if not (0.0 < length_fraction <= 1.0):
        raise DesignError("length_fraction must be in (0, 1]")
    X = result.X
    L = float(X[-1] - X[0])
    x_cut = X[0] + length_fraction * L
    idx = int(np.searchsorted(X, x_cut, side="right") - 1)
    idx = max(0, min(idx, len(X) - 1))

    Isp_t = float(result.Isp_cum[idx])
    F_t = Isp_t * result.m_dot * G0
    Cf_t = F_t / (result.Pc * result.At)

    sl = slice(0, idx + 1)
    contour = Contour(
        x=result.X[sl].copy(),
        r=result.R[sl].copy(),
        M=result.M[sl].copy(),
    )
    return SpikeResult(
        gamma=result.gamma,
        molar_mass=result.molar_mass,
        Pc=result.Pc,
        Tc=result.Tc,
        Pa=result.Pa,
        Re=result.Re,
        expansion_ratio=result.expansion_ratio,
        exit_mach=result.exit_mach,
        delta_rad=result.delta_rad,
        ht=result.ht,
        At=result.At,
        m_dot=result.m_dot,
        F=F_t,
        Isp=Isp_t,
        Cf=Cf_t,
        Pe=result.Pe,
        X=result.X[sl].copy(),
        R=result.R[sl].copy(),
        M=result.M[sl].copy(),
        P=result.P[sl].copy(),
        T=result.T[sl].copy(),
        Isp_cum=result.Isp_cum[sl].copy(),
        contour=contour,
        meta={
            **result.meta,
            "truncated": True,
            "length_fraction": length_fraction,
            "full_Isp_s": result.Isp,
        },
    )
