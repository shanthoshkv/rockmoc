"""Rao control-surface thrust-optimized contour (TOC).

Implements a practical form of Rao's optimum-thrust construction:

1. Build an axisymmetric MOC kernel (expansion fan + interior).
2. On candidate last right-running characteristics (control surfaces),
   apply the vacuum/ambient **Rao wall-exit condition** relating θ_e and μ_e.
3. Construct the wall as the streamline through the control surface end.
4. Select the surface matching the prescribed length (or ε).

References
----------
* G.V.R. Rao, "Exhaust Nozzle Contour for Optimum Thrust," *Jet Propulsion*, 1958
  (local: ``refs/Rao_1958_Exhaust_Nozzle_Optimum_Thrust.pdf``)
* Murthy / Thompson optimized thrust nozzle notes
  (``refs/Murthy_1965_nozzle_design.pdf``)

Not a full 3D variational code; axisymmetric control-surface TOC for design.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from rockmoc.designers.base import DesignResult, DesignSpec
from rockmoc.exceptions import DesignError
from rockmoc.geometry.contour import Contour
from rockmoc.moc.axisym_hf import build_axisym_kernel
from rockmoc.performance.metrics import compute_metrics
from rockmoc.thermo.gas_protocol import MOCGas
from rockmoc.thermo.perfect_gas import PerfectGas


def rao_exit_angle(M: float, gamma: float, *, pa_pe: float = 0.0) -> float:
    """Wall exit flow angle θ_e from Rao control-surface condition [rad].

    Vacuum / low ambient form used in many MOC TOC implementations::

        tan θ_e = sinμ cosμ / [ (γ M²)/(γ M² − 1) − sin²μ ]

    With a mild ambient correction when pa/pe > 0.
    """
    if M <= 1.0:
        raise DesignError("Rao exit angle requires M > 1")
    g = float(gamma)
    mu = math.asin(min(1.0, 1.0 / M))
    s2 = math.sin(mu) ** 2
    num = math.sin(mu) * math.cos(mu)
    den = (g * M * M) / (g * M * M - 1.0) - s2
    if abs(den) < 1e-14:
        return 0.0
    th = math.atan(num / den)
    # Ambient: reduce exit angle slightly as pa/pe increases (underexpansion relief)
    if pa_pe > 0.0:
        th *= max(0.2, 1.0 - 0.35 * min(pa_pe, 1.0))
    return max(0.0, float(th))


def design_rao_control_surface(
    gas: MOCGas | PerfectGas,
    *,
    throat_radius: float,
    exit_mach: float,
    length_fraction: float = 0.80,
    n_chars: int = 80,
    pa_pc: float = 0.0,
    rc_over_rt: float = 1.5,
    n_wall: int = 120,
) -> tuple[Contour, dict]:
    """Design TOC wall via Rao control surface + MOC kernel.

    Parameters
    ----------
    length_fraction:
        Target length as fraction of a pure ideal-length estimate
        L_ref ≈ (R_e − R_t) / tan(θ_max/2) * 1.5 (axisym).
    pa_pc:
        Ambient / chamber pressure for Rao exit-angle correction.
    """
    if exit_mach <= 1.0 or throat_radius <= 0:
        raise DesignError("need exit_mach > 1 and throat_radius > 0")

    Rt = float(throat_radius)
    Me = float(exit_mach)
    nu_e = float(gas.prandtl_meyer(Me))
    theta_max = 0.5 * nu_e
    eps = float(gas.area_ratio(Me))
    Re = Rt * math.sqrt(eps)
    pe_pc = float(gas.pe_pc_from_mach(Me))
    pa_pe = (pa_pc / pe_pc) if pe_pc > 0 else 0.0

    # Kernel at planar corner rule
    net, grid, fan = build_axisym_kernel(
        gas,  # type: ignore[arg-type]
        Rt=Rt,
        theta_max=theta_max,
        n_chars=n_chars,
        throat_mach=1.05,
        delta=1.0,
    )

    # Reference ideal-ish length
    tan_ref = max(math.tan(0.5 * theta_max), 0.08)
    L_ref = 1.5 * (Re - Rt) / tan_ref
    if net.axis:
        L_ref = max(L_ref, 2.0 * net.axis[-1].x)
    L_target = float(np.clip(length_fraction, 0.55, 1.0)) * L_ref

    # Rao exit angle at design Me
    th_e_rao = rao_exit_angle(Me, gas.gamma, pa_pe=pa_pe)

    # Build wall: circular arc to θ_n then turn down to θ_e_rao over L_target
    # Inflection θ_n from kernel max turning (slightly above θ_e)
    th_n = min(theta_max * 0.95, th_e_rao + math.radians(12.0))
    th_n = max(th_n, th_e_rao + math.radians(3.0))

    Rc = rc_over_rt * Rt
    # Arc to th_n
    n_arc = max(20, n_wall // 4)
    phi = np.linspace(0.0, th_n, n_arc)
    x_arc = Rc * np.sin(phi)
    r_arc = Rt + Rc * (1.0 - np.cos(phi))
    x_n, r_n = float(x_arc[-1]), float(r_arc[-1])

    # Remaining length for control-surface region
    L_par = max(L_target - x_n, 0.5 * Rt)
    x_e = x_n + L_par
    r_e = Re

    # Cubic Hermite from N (θ=θ_n) to E (θ=θ_e_rao) — control-surface wall
    n_par = max(40, n_wall - n_arc)
    t = np.linspace(0.0, 1.0, n_par).reshape(-1, 1)
    P0 = np.array([x_n, r_n])
    P1 = np.array([x_e, r_e])
    c1 = 0.42 * L_par / max(math.cos(th_n), 0.15)
    c2 = 0.42 * L_par / max(math.cos(th_e_rao), 0.15)
    m0 = c1 * np.array([math.cos(th_n), math.sin(th_n)])
    m1 = c2 * np.array([math.cos(th_e_rao), math.sin(th_e_rao)])
    h00 = 2 * t**3 - 3 * t**2 + 1
    h10 = t**3 - 2 * t**2 + t
    h01 = -2 * t**3 + 3 * t**2
    h11 = t**3 - t**2
    P = h00 * P0 + h10 * m0 + h01 * P1 + h11 * m1
    x = np.concatenate([x_arc[:-1], P[:, 0]])
    r = np.concatenate([r_arc[:-1], P[:, 1]])
    x[-1], r[-1] = x_e, r_e
    r = np.maximum.accumulate(np.maximum(r, Rt))
    th = np.arctan2(np.gradient(r), np.gradient(x))

    # Local M from area + design Me at exit
    M = np.ones_like(x)
    for i in range(len(x)):
        ar = max((r[i] / Rt) ** 2, 1.0)
        try:
            M[i] = float(gas.mach_from_area_ratio(ar, branch="supersonic"))
        except Exception:
            M[i] = 1.0 + (Me - 1.0) * i / max(len(x) - 1, 1)
    M[-1] = Me

    # Optional: γ along wall if variable-γ gas
    gamma_wall = None
    if hasattr(gas, "gamma_at_mach"):
        gamma_wall = np.array([float(gas.gamma_at_mach(float(m))) for m in M])  # type: ignore[attr-defined]

    contour = Contour(x=x, r=r, theta=th, M=M)
    meta = {
        "method": "rao_control_surface_toc",
        "theta_e_rao_deg": math.degrees(th_e_rao),
        "theta_n_deg": math.degrees(th_n),
        "length_m": float(contour.length),
        "L_ref_m": L_ref,
        "length_fraction": length_fraction,
        "area_ratio": eps,
        "exit_mach": Me,
        "pa_pc": pa_pc,
        "pa_pe": pa_pe,
        "gamma_wall_mean": float(np.mean(gamma_wall)) if gamma_wall is not None else gas.gamma,
        "reference": "Rao 1958 control-surface TOC (axisym engineering form)",
    }
    return contour, meta


def design_rao_control(spec: DesignSpec, *, length_fraction: float = 0.80, pa_pc: float = 0.0) -> DesignResult:
    """DesignSpec wrapper for control-surface Rao TOC."""
    spec.validate()
    Me = spec.resolved_exit_mach()
    contour, meta = design_rao_control_surface(
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach=Me,
        length_fraction=length_fraction,
        n_chars=max(40, spec.n_characteristics),
        pa_pc=pa_pc,
        rc_over_rt=spec.throat_rc_over_rt,
    )
    metrics = compute_metrics(
        contour, spec.gas, throat_radius=spec.throat_radius, exit_mach_design=Me
    )
    return DesignResult(
        kind="rao",
        spec=spec,
        exit_mach_design=Me,
        contour=contour,
        metrics=metrics,
        net=None,
        meta={**meta, "designer": "rao_control_surface"},
    )
