"""Free-shock flow separation in over-expanded nozzles.

Criteria give the incipient-separation pressure ratio p_sep/p_a (all from the
cited sources; scatter is large, NASA advised a 20 % margin):

=====================  ==========================================  =====================
``summerfield``        0.4                                         conical, 1950s
``kalt_badal``         (2/3)(p0/pa)^−0.2                           RocketCEA default
``schilling``          0.582 (p0/pa)^−0.195 (contoured)            Östlund 2004
                       0.541 (p0/pa)^−0.136 (conical)
``schmucker``          (1.88 M − 1)^−0.64                          European default
``stark``              π / (3 M)                                   Stark, AIAA 2005-3940
``morrisette``         1.082 − 0.363 M + 0.0386 M², 2.4 ≤ M ≤ 4.5  NASA TP-1207
=====================  ==========================================  =====================

Separated thrust: the nozzle acts as if truncated at ε_sep (``recovery="none"``,
wall pressure ≈ p_a downstream) or, with ``recovery="kalt_badal"``, with the
Kalt & Badal pressure-recovery integral exactly as RocketCEA ``separated_Cf``.
NASA TP-1207 guards are reported: full flow when p0/pa ≥ 1 + 0.39 p0/pe, and
criteria lose validity in the end-effect zone ε_s > 0.8(ε_e − 1) + 1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from rockmoc.exceptions import DesignError
from rockmoc.thermo.gas_protocol import MOCGas

Criterion = Literal["summerfield", "kalt_badal", "schilling", "schmucker", "stark", "morrisette"]


def psep_over_pa(criterion: Criterion, *, M: float, p0_pa: float, conical: bool = False) -> float:
    """Incipient-separation pressure / ambient for wall Mach ``M`` and p0/pa."""
    if criterion == "summerfield":
        return 0.4
    if criterion == "kalt_badal":
        return (2.0 / 3.0) * p0_pa**-0.2
    if criterion == "schilling":
        return 0.541 * p0_pa**-0.136 if conical else 0.582 * p0_pa**-0.195
    if criterion == "schmucker":
        return (1.88 * M - 1.0) ** -0.64
    if criterion == "stark":
        return math.pi / (3.0 * M)
    if criterion == "morrisette":
        return 1.082 - 0.363 * M + 0.0386 * M * M
    raise DesignError(f"unknown separation criterion {criterion!r}")


@dataclass(slots=True)
class SeparationResult:
    separated: bool
    eps_sep: float  # area ratio at separation (= ε_e when attached)
    M_sep: float
    psep_pa: float
    Cf: float  # ambient thrust coefficient (1-D)
    Cf_attached: float  # what attached-flow theory would give
    full_flow_tp1207: bool
    end_effect: bool
    criterion: str


def separated_cf(
    gas: MOCGas,
    *,
    area_ratio: float,
    pc: float,
    pa: float,
    criterion: Criterion = "schmucker",
    recovery: Literal["none", "kalt_badal"] = "none",
    conical: bool = False,
) -> SeparationResult:
    """1-D thrust coefficient of an over-expanded nozzle with free-shock separation."""
    if pa <= 0.0:
        Me = float(gas.mach_from_area_ratio(area_ratio, branch="supersonic"))
        cf = float(gas.thrust_coefficient(float(gas.pe_pc_from_mach(Me)), 0.0, exit_mach=Me, area_ratio=area_ratio))
        return SeparationResult(False, area_ratio, Me, 0.0, cf, cf, True, False, criterion)
    p0_pa = pc / pa

    def M_of(eps: float) -> float:
        return float(gas.mach_from_area_ratio(eps, branch="supersonic"))

    def cf_at(eps: float, pa_pc: float) -> float:
        M = M_of(eps)
        return float(gas.thrust_coefficient(float(gas.pe_pc_from_mach(M)), pa_pc, exit_mach=M, area_ratio=eps))

    def resid(eps: float) -> float:  # wall p/pa minus criterion (positive while attached)
        M = M_of(eps)
        return float(gas.pe_pc_from_mach(M)) * p0_pa - psep_over_pa(criterion, M=M, p0_pa=p0_pa, conical=conical)

    Me = M_of(area_ratio)
    cf_att = cf_at(area_ratio, pa / pc)
    pe_pc = float(gas.pe_pc_from_mach(Me))
    full = p0_pa >= 1.0 + 0.39 / pe_pc
    if resid(area_ratio) >= 0.0:
        return SeparationResult(False, area_ratio, Me, psep_over_pa(criterion, M=Me, p0_pa=p0_pa, conical=conical),
                                cf_att, cf_att, full, False, criterion)
    lo, hi = 1.0 + 1e-9, area_ratio
    if resid(lo) <= 0.0:
        raise DesignError("flow separates at the throat: p0/pa too low for a supersonic nozzle")
    for _ in range(100):  # bisection on ε (resid monotone in practice)
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if resid(mid) > 0.0 else (lo, mid)
    eps_s = 0.5 * (lo + hi)
    Ms = M_of(eps_s)
    ps = float(gas.pe_pc_from_mach(Ms)) * pc
    if recovery == "kalt_badal":  # RocketCEA separated_Cf.sepNozzleCf
        c = 2.4 / 1.45
        if eps_s <= area_ratio * c / (1 + c) + 1 / (1 + c):
            e95 = eps_s + (eps_s - 1) / 2.4
        else:
            e95 = eps_s + (area_ratio - eps_s)
        cf = (cf_at(eps_s, 0.0) + 0.55 * (ps + 0.95 * pa) * (e95 - eps_s) / pc
              - (0.975 * e95 + 0.025 * area_ratio) * pa / pc)
    else:
        cf = cf_at(eps_s, pa / pc)
    return SeparationResult(True, eps_s, Ms, ps / pa, cf, cf_att, full,
                            eps_s > 0.8 * (area_ratio - 1.0) + 1.0, criterion)
