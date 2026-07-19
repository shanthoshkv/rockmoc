"""Design sweeps for performance charts (OF, gamma, area ratio)."""

from __future__ import annotations

from typing import Any

import numpy as np

from rockmoc.designers.base import DesignSpec
from rockmoc.designers.mln import design_mln
from rockmoc.performance.integrated import integrated_performance, off_design_altitude_sweep
from rockmoc.thermo.perfect_gas import PerfectGas


def altitude_sweep_rows(
    contour: Any,
    gas: PerfectGas,
    *,
    throat_radius: float,
    exit_mach: float,
    pc: float = 5e6,
    Tc: float = 3200.0,
    altitudes_m: list[float] | None = None,
) -> list[dict]:
    return off_design_altitude_sweep(
        contour,
        gas,
        throat_radius=throat_radius,
        exit_mach=exit_mach,
        pc=pc,
        Tc=Tc,
        altitudes_m=altitudes_m,
    )


def gamma_sweep_rows(
    *,
    throat_radius: float = 0.05,
    exit_mach: float = 3.0,
    gammas: list[float] | None = None,
    n_chars: int = 24,
    R: float = 350.0,
    pc: float = 5e6,
    Tc: float = 3200.0,
) -> list[dict]:
    if gammas is None:
        gammas = list(np.linspace(1.15, 1.35, 9))
    rows = []
    for g in gammas:
        gas = PerfectGas(float(g), R)
        res = design_mln(
            DesignSpec(
                gas=gas,
                throat_radius=throat_radius,
                exit_mach=exit_mach,
                n_characteristics=n_chars,
            ),
            store_net=False,
        )
        perf = integrated_performance(
            res.contour, gas, throat_radius=throat_radius, exit_mach=exit_mach, pc=pc, Tc=Tc, pa=0.0
        )
        rows.append(
            {
                "gamma": float(g),
                "L_over_Rt": res.metrics.length / throat_radius,
                "area_ratio": res.metrics.area_ratio,
                "Cf_vac": perf.Cf_vac,
                "Isp_vac_s": perf.Isp_vac,
            }
        )
    return rows


def of_ratio_sweep_rows(
    *,
    throat_radius: float = 0.05,
    exit_mach: float = 3.0,
    of_list: list[float] | None = None,
    n_chars: int = 20,
    pc: float = 5e6,
    Tc: float = 3200.0,
    use_cea: bool = True,
    ox: str = "LOX",
    fuel: str = "RP1",
    Pc_psia: float = 725.0,
) -> list[dict]:
    """Sweep oxidizer/fuel ratio.

    If RocketCEA is available and ``use_cea``, use real chamber γ and Tc.
    Otherwise map OF → γ with an engineering correlation for demo plots.
    """
    if of_list is None:
        of_list = list(np.linspace(1.8, 3.2, 8))

    cea = None
    if use_cea:
        try:
            from rocketcea.cea_obj import CEA_Obj  # type: ignore

            cea = CEA_Obj(oxName=ox, fuelName=fuel)
        except Exception:
            cea = None

    rows = []
    for of in of_list:
        of = float(of)
        if cea is not None:
            try:
                mw, gamma = cea.get_Chamber_MolWt_gamma(Pc_psia, of, 40.0)
                try:
                    _, gamma_t = cea.get_Throat_MolWt_gamma(Pc_psia, of, 40.0)
                    gamma = gamma_t
                except Exception:
                    pass
                Tc_use = float(cea.get_Tcomb(Pc_psia, of))
                R = 8314.462618 / float(mw)
                g = float(gamma)
            except Exception:
                g, R, Tc_use = _synthetic_of_props(of)
        else:
            g, R, Tc_use = _synthetic_of_props(of)
            Tc_use = Tc if Tc_use is None else Tc_use

        gas = PerfectGas(g, R)
        res = design_mln(
            DesignSpec(
                gas=gas,
                throat_radius=throat_radius,
                exit_mach=exit_mach,
                n_characteristics=n_chars,
            ),
            store_net=False,
        )
        perf = integrated_performance(
            res.contour,
            gas,
            throat_radius=throat_radius,
            exit_mach=exit_mach,
            pc=pc,
            Tc=Tc_use or Tc,
            pa=0.0,
        )
        rows.append(
            {
                "OF": of,
                "gamma": g,
                "Tc_K": Tc_use or Tc,
                "area_ratio": res.metrics.area_ratio,
                "L_over_Rt": res.metrics.length / throat_radius,
                "Cf_vac": perf.Cf_vac,
                "Isp_vac_s": perf.Isp_vac,
                "cea": cea is not None,
            }
        )
    return rows


def _synthetic_of_props(of: float) -> tuple[float, float, float]:
    """Engineering stand-in for LOX/RP-1 style OF → (γ, R, Tc)."""
    # Stoichiometric ~2.3; γ min near stoich, Tc peaks near stoich
    of0 = 2.3
    g = 1.22 + 0.04 * abs(of - of0)  # higher γ off-stoich
    g = float(np.clip(g, 1.15, 1.35))
    R = 360.0 - 15.0 * abs(of - of0)
    Tc = 3400.0 - 400.0 * (of - of0) ** 2
    return g, R, Tc


def area_ratio_sweep_rows(
    *,
    throat_radius: float = 0.05,
    gamma: float = 1.25,
    eps_list: list[float] | None = None,
    n_chars: int = 20,
    pc: float = 5e6,
    Tc: float = 3200.0,
) -> list[dict]:
    if eps_list is None:
        eps_list = list(np.geomspace(4.0, 60.0, 10))
    gas = PerfectGas(gamma)
    rows = []
    for eps in eps_list:
        res = design_mln(
            DesignSpec(
                gas=gas,
                throat_radius=throat_radius,
                area_ratio=float(eps),
                n_characteristics=n_chars,
            ),
            store_net=False,
        )
        Me = res.exit_mach_design
        perf = integrated_performance(
            res.contour, gas, throat_radius=throat_radius, exit_mach=Me, pc=pc, Tc=Tc, pa=0.0
        )
        rows.append(
            {
                "area_ratio_target": float(eps),
                "area_ratio": res.metrics.area_ratio,
                "exit_mach": Me,
                "L_over_Rt": res.metrics.length / throat_radius,
                "Cf_vac": perf.Cf_vac,
                "Isp_vac_s": perf.Isp_vac,
            }
        )
    return rows
