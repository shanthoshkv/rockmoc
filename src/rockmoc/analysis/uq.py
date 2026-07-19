"""Uncertainty quantification for nozzle design metrics.

Monte Carlo / one-at-a-time sensitivity on γ, mesh density, and optional BL.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np

from rockmoc.designers.base import DesignSpec
from rockmoc.designers.mln import design_mln
from rockmoc.exceptions import DesignError
from rockmoc.thermo.perfect_gas import PerfectGas


@dataclass
class UQResult:
    """Summary statistics for sampled designs."""

    samples: list[dict[str, Any]]
    stats: dict[str, dict[str, float]]
    meta: dict

    def summary(self) -> dict:
        return {"n_samples": len(self.samples), "stats": self.stats, **self.meta}


def _stats(arr: np.ndarray) -> dict[str, float]:
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "p05": float(np.percentile(arr, 5)),
        "p95": float(np.percentile(arr, 95)),
    }


def monte_carlo_mln(
    *,
    throat_radius: float,
    exit_mach: float,
    gamma_mean: float = 1.25,
    gamma_std: float = 0.02,
    n_chars_choices: list[int] | None = None,
    n_samples: int = 24,
    seed: int = 42,
    R: float = 350.0,
) -> UQResult:
    """Sample γ and mesh density; collect L/Rt, ε, Cf_vac (1-D)."""
    if n_samples < 3:
        raise DesignError("n_samples must be ≥ 3")
    rng = np.random.default_rng(seed)
    if n_chars_choices is None:
        n_chars_choices = [24, 40, 60]

    samples: list[dict[str, Any]] = []
    for i in range(n_samples):
        g = float(rng.normal(gamma_mean, gamma_std))
        g = float(np.clip(g, 1.08, 1.65))
        n = int(rng.choice(n_chars_choices))
        gas = PerfectGas(gamma=g, R=R)
        try:
            res = design_mln(
                DesignSpec(
                    gas=gas,
                    throat_radius=throat_radius,
                    exit_mach=exit_mach,
                    n_characteristics=n,
                ),
                store_net=False,
            )
            pe_pc = float(gas.pe_pc_from_mach(exit_mach))
            eps = res.metrics.area_ratio
            Cf = gas.thrust_coefficient(pe_pc, 0.0, exit_mach=exit_mach, area_ratio=eps)
            samples.append(
                {
                    "gamma": g,
                    "n_chars": n,
                    "L_over_Rt": res.metrics.length / throat_radius,
                    "area_ratio": eps,
                    "Cf_vac_1d": float(Cf),
                    "ok": True,
                }
            )
        except Exception as exc:
            samples.append(
                {
                    "gamma": g,
                    "n_chars": n,
                    "ok": False,
                    "error": str(exc),
                }
            )

    ok = [s for s in samples if s.get("ok")]
    if not ok:
        raise DesignError("All UQ samples failed")
    stats = {
        "L_over_Rt": _stats(np.array([s["L_over_Rt"] for s in ok])),
        "area_ratio": _stats(np.array([s["area_ratio"] for s in ok])),
        "Cf_vac_1d": _stats(np.array([s["Cf_vac_1d"] for s in ok])),
        "gamma": _stats(np.array([s["gamma"] for s in ok])),
    }
    return UQResult(
        samples=samples,
        stats=stats,
        meta={
            "n_ok": len(ok),
            "n_fail": len(samples) - len(ok),
            "exit_mach": exit_mach,
            "throat_radius": throat_radius,
            "gamma_mean": gamma_mean,
            "gamma_std": gamma_std,
            "seed": seed,
        },
    )


def sensitivity_one_at_a_time(
    *,
    throat_radius: float,
    exit_mach: float,
    gamma: float = 1.25,
    R: float = 350.0,
    n_chars: int = 40,
    d_gamma: float = 0.02,
) -> dict[str, Any]:
    """Finite-difference sensitivity of L/Rt and ε to γ."""
    def run(g: float) -> dict[str, float]:
        res = design_mln(
            DesignSpec(
                gas=PerfectGas(g, R),
                throat_radius=throat_radius,
                exit_mach=exit_mach,
                n_characteristics=n_chars,
            ),
            store_net=False,
        )
        return {
            "L_over_Rt": res.metrics.length / throat_radius,
            "area_ratio": res.metrics.area_ratio,
        }

    base = run(gamma)
    hi = run(gamma + d_gamma)
    lo = run(max(gamma - d_gamma, 1.1))
    return {
        "base": base,
        "d_gamma": d_gamma,
        "dL_dgamma": (hi["L_over_Rt"] - lo["L_over_Rt"]) / (2 * d_gamma),
        "deps_dgamma": (hi["area_ratio"] - lo["area_ratio"]) / (2 * d_gamma),
        "hi": hi,
        "lo": lo,
    }


def mesh_certificate(
    *,
    throat_radius: float,
    exit_mach: float,
    gamma: float = 1.25,
    n_coarse: int = 24,
    n_fine: int = 48,
    tol_rel: float = 0.02,
) -> dict[str, Any]:
    """Check relative change in L and ε between coarse and fine mesh."""
    gas = PerfectGas(gamma)
    r0 = design_mln(
        DesignSpec(
            gas=gas,
            throat_radius=throat_radius,
            exit_mach=exit_mach,
            n_characteristics=n_coarse,
        ),
        store_net=False,
    )
    r1 = design_mln(
        DesignSpec(
            gas=gas,
            throat_radius=throat_radius,
            exit_mach=exit_mach,
            n_characteristics=n_fine,
        ),
        store_net=False,
    )
    dL = abs(r1.metrics.length - r0.metrics.length) / max(r0.metrics.length, 1e-12)
    deps = abs(r1.metrics.area_ratio - r0.metrics.area_ratio) / max(
        r0.metrics.area_ratio, 1e-12
    )
    passed = dL < tol_rel and deps < tol_rel
    return {
        "passed": passed,
        "tol_rel": tol_rel,
        "dL_rel": dL,
        "deps_rel": deps,
        "n_coarse": n_coarse,
        "n_fine": n_fine,
        "L_coarse": r0.metrics.length,
        "L_fine": r1.metrics.length,
        "eps_coarse": r0.metrics.area_ratio,
        "eps_fine": r1.metrics.area_ratio,
    }
