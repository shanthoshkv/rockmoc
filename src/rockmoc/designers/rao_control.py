"""Rao thrust-optimised contour (TOC) — DesignSpec wrapper.

The contour is built by :func:`rockmoc.moc.rao.design_rao_moc`: smooth-throat
characteristic lattice to the inflection N, Rao's control surface DE
(V cos(θ−μ)/cos μ and r ρV² sin²θ tan μ constant, mass-matched), and the wall
NE as the mass-flow streamline between ND and DE.

References
----------
* G.V.R. Rao, "Exhaust Nozzle Contour for Optimum Thrust," *Jet Propulsion*, 1958
* Zucrow & Hoffman, *Gas Dynamics* Vol. 2, §16.6
"""

from __future__ import annotations

import math

from rockmoc.designers.base import DesignResult, DesignSpec
from rockmoc.exceptions import DesignError
from rockmoc.geometry.contour import Contour
from rockmoc.moc.rao import design_rao_moc
from rockmoc.moc.types import CharNet
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
    rc_over_rt: float = 0.382,
    throat_mach: float = 1.02,
) -> tuple[Contour, dict, CharNet]:
    """Rao TOC by MOC (:func:`rockmoc.moc.rao.design_rao_moc`), vacuum optimum.

    ``length_fraction`` is relative to the 15° cone (Rao's convention). An
    unreachable length falls back to the nearest feasible one with a warning.
    ``pa_pc`` > 0 is reported only: the construction is the vacuum optimum.
    N ≥ 64 and a near-sonic start (M₀ = 1.02) keep the short-bell end of the
    feasible length range resolved.
    """
    contour, net, meta = design_rao_moc(
        gas, throat_radius=throat_radius, exit_mach=exit_mach,
        length_fraction=length_fraction, n_chars=max(int(n_chars), 64),
        rc_over_rt=rc_over_rt, throat_mach=throat_mach, clamp_length=True,
    )
    pe_pc = float(gas.pe_pc_from_mach(exit_mach))
    pa_pe = (pa_pc / pe_pc) if pe_pc > 0 else 0.0
    if pa_pc > 0.0:
        meta["warnings"].append("pa_pc > 0 ignored: contour is the vacuum Rao optimum")  # input note, not a fallback
    meta.update(
        theta_e_rao_deg=math.degrees(rao_exit_angle(exit_mach, gas.gamma, pa_pe=pa_pe)),
        pa_pc=pa_pc,
        pa_pe=pa_pe,
        exit_mach=exit_mach,
    )
    return contour, meta, net


def design_rao_control(spec: DesignSpec, *, length_fraction: float = 0.80, pa_pc: float = 0.0) -> DesignResult:
    """DesignSpec wrapper for control-surface Rao TOC."""
    spec.validate()
    Me = spec.resolved_exit_mach()
    contour, meta, net = design_rao_control_surface(
        spec.gas,
        throat_radius=spec.throat_radius,
        exit_mach=Me,
        length_fraction=length_fraction,
        n_chars=max(40, spec.n_characteristics),
        pa_pc=pa_pc,
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
        net=net,
        meta={**meta, "designer": "rao_control_surface"},
    )
