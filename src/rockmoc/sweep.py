"""Batch parameter sweeps: Cartesian grid over any design/operating parameter.

Axes are ``name -> list of values``; parse CLI-style strings with :func:`parse_axis`::

    "Me=2:4:5"            # start:stop:count  (linspace)
    "area_ratio=10,25,40" # explicit list
    "type=mln,rao"        # categorical
    "pc_bar=30:100:8"     # unit suffix _bar/_mpa/_kpa/_psia on pressures

Parameters
----------
design : ``Me`` ``area_ratio`` ``pe_pc`` (pe/pc) ``pe`` (exit static, Pa; needs ``pc``)
         ``Rt`` (throat radius, in the base spec's ``[design] units``) ``n`` ``rc`` ``type``
         ``length`` (Rao) ``truncate`` (TIC)
gas      : ``gamma`` ``R`` (constant-γ gas) ``of`` (mixture ratio; needs ``[propellants]``)
           ``gas_model`` = perfect | equilibrium | frozen (real-gas RocketCEA tables;
           sweep ``gas_model=frozen,equilibrium`` to bracket real engine performance)
operating: ``pc`` ``Tc`` ``pa`` ``altitude`` (m)

Each case yields one row (geometry + Cf/Isp/thrust/mdot); a failing case records
``error`` and the sweep continues.
"""

from __future__ import annotations

import copy
import csv
import itertools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

from rockmoc.exceptions import DesignError

_PRESSURE_UNITS = {"_bar": 1e5, "_mpa": 1e6, "_kpa": 1e3, "_psia": 6894.757293168}
_EXIT_TARGETS = ("Me", "area_ratio", "pe_pc", "pe")
PARAMS = (*_EXIT_TARGETS, "Rt", "n", "rc", "type", "length", "truncate",
          "gamma", "R", "of", "gas_model", "pc", "Tc", "pa", "altitude")


def parse_axis(text: str) -> tuple[str, list]:
    """``"Me=2:4:5"`` / ``"Me=2,3,4"`` / ``"type=mln,rao"`` → (name, values)."""
    name, _, rhs = text.partition("=")
    name = name.strip()
    if not rhs:
        raise DesignError(f"bad --vary {text!r}; use NAME=start:stop:count or NAME=a,b,c")
    if ":" in rhs:
        a, b, n = rhs.split(":")
        vals: list = [float(v) for v in np.linspace(float(a), float(b), int(n))]
    else:
        vals = []
        for v in rhs.split(","):
            try:
                vals.append(float(v))
            except ValueError:
                vals.append(v.strip())
    return name, vals


def _canon(name: str, vals: list) -> tuple[str, list]:
    """Strip a pressure-unit suffix (``pc_bar`` → ``pc`` in Pa); validate the name."""
    for suf, k in _PRESSURE_UNITS.items():
        if name.endswith(suf):
            return name[: -len(suf)], [v * k for v in vals]
    return name, vals


def _apply(spec: dict, case: dict) -> dict:
    """Return a copy of ``spec`` with one sweep case applied."""
    s = copy.deepcopy(spec)
    d = s.setdefault("design", {})
    perf = s.setdefault("performance", {})
    g = s.setdefault("gas", {})
    for k, v in case.items():
        if k in ("Me", "area_ratio", "pe_pc"):
            for t in ("exit_mach", "area_ratio", "pressure_ratio"):
                d.pop(t, None)
            d[{"Me": "exit_mach", "area_ratio": "area_ratio", "pe_pc": "pressure_ratio"}[k]] = v
        elif k == "pe":
            pass  # resolved below (needs final pc)
        elif k == "Rt":
            d["throat_radius"] = v
        elif k == "n":
            d["n_characteristics"] = int(v)
        elif k == "rc":
            d["throat_rc_over_rt"] = v
        elif k == "type":
            d["type"] = v
        elif k == "length":
            d["length_fraction"] = v
        elif k == "truncate":
            d["truncate_fraction"] = v
        elif k in ("gamma", "R"):
            g[k] = v
            s.pop("propellants", None)  # explicit gas beats propellants
        elif k == "gas_model":
            g["model"] = v
        elif k == "of":
            if "propellants" not in s:
                raise DesignError("sweeping 'of' needs a [propellants] ox/fuel in the base spec")
            s["propellants"]["of"] = v
        elif k in ("pc", "Tc", "pa"):
            perf[k] = v
        elif k == "altitude":
            perf["altitude"] = v
            perf.pop("pa", None)
        else:
            raise DesignError(f"unknown sweep parameter {k!r}; choose from {PARAMS}")
    if "propellants" in s:
        s["propellants"]["pc"] = perf.get("pc", s["propellants"].get("pc", 5e6))
    if "pe" in case:
        for t in ("exit_mach", "area_ratio"):
            d.pop(t, None)
        d["pressure_ratio"] = case["pe"] / perf.get("pc", 5e6)
    return s


def _run_case(args: tuple) -> dict:
    """Worker (top-level → picklable for process pools)."""
    idx, spec, case, export_dir, export_formats = args
    from rockmoc.spec_io import (
        build_design,
        export_design,
        performance_from_spec,
        validate_spec,
    )

    row: dict[str, Any] = {"case": idx, **case}
    try:
        validate_spec(spec)
        res = build_design(spec)
        perf = performance_from_spec(res, spec)  # _apply guarantees [performance]
        m = res.metrics
        Rt = res.spec.throat_radius
        gas = res.spec.gas
        row.update(
            kind=res.kind, area_ratio=m.area_ratio, Me=res.exit_mach_design, pe_pc=m.pe_pc_1d,
            gas_model=getattr(getattr(gas, "table", None), "mode", "perfect"),
            gamma_throat=gas.gamma_at_mach(1.0),
            gamma_exit=gas.gamma_at_mach(float(res.exit_mach_design)),
            R_used=gas.R, Rt_m=Rt, Re_m=m.exit_radius, L_m=m.length,
            L_over_Rt=m.length / Rt, theta_max_deg=float(np.degrees(m.max_wall_angle)),
            theta_exit_deg=perf["theta_exit_deg"], pc_Pa=perf["pc_Pa"], Tc_K=perf["Tc_K"],
            pa_Pa=perf["pa_Pa"], cstar_m_s=perf["cstar_m_s"], Cf_vac=perf["Cf_vac"],
            Cf_amb=perf["Cf_amb"], Isp_vac_s=perf["Isp_vac_s"], Isp_amb_s=perf["Isp_amb_s"],
            thrust_vac_N=perf["thrust_vac_N"], thrust_amb_N=perf["thrust_amb_N"],
            mdot_kg_s=perf["mdot_kg_s"], eta_nozzle=perf["nozzle_efficiency"],
        )
        if "propellants" in spec and row["gas_model"] == "perfect":
            row["cea"] = _used_cea(spec)
        if export_dir and export_formats:
            e = {**spec.get("export", {}), "formats": export_formats, "out_dir": str(export_dir),
                 "name": f"case_{idx:04d}"}
            export_design(res, {**spec, "export": e})
    except Exception as exc:  # noqa: BLE001  # ponytail: one bad case must not kill the grid
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row


def _used_cea(spec: dict) -> bool:
    from rockmoc.thermo.propellant import of_props

    pr = spec["propellants"]
    return of_props(pr.get("ox", "LOX"), pr.get("fuel", "RP1"), pr["of"],
                    pr.get("pc", 5e6))["cea"]


def run_sweep(
    base: dict,
    axes: dict[str, list],
    *,
    export_dir: str | Path | None = None,
    export_formats: list[str] | None = None,
    jobs: int = 1,
) -> list[dict]:
    """Run the Cartesian product of ``axes`` over the ``base`` spec dict (TOML layout).

    Returns one row dict per case. ``export_formats`` (e.g. ``["csv"]``) also writes each
    case's wall to ``export_dir/case_NNNN.*`` (units/points from ``base["export"]``).
    """
    canon = dict(_canon(k, v) for k, v in axes.items())
    hit = [k for k in canon if k in _EXIT_TARGETS]
    if len(hit) > 1:
        raise DesignError(f"sweep only one exit target at a time; got {hit}")
    for k in canon:
        if k not in PARAMS:
            raise DesignError(f"unknown sweep parameter {k!r}; choose from {PARAMS}")
    if not canon:
        raise DesignError("no sweep axes given (use --vary NAME=start:stop:count)")
    names = list(canon)
    cases = [dict(zip(names, vals, strict=True)) for vals in itertools.product(*canon.values())]
    if export_dir:
        Path(export_dir).mkdir(parents=True, exist_ok=True)
    work = [(i, _apply(base, c), c, export_dir, export_formats) for i, c in enumerate(cases, 1)]
    if jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            return list(ex.map(_run_case, work))
    return [_run_case(w) for w in work]


def write_rows(rows: list[dict], path: str | Path | None = None) -> str:
    """CSV of sweep rows; writes to ``path`` if given. Returns the CSV text."""
    import io

    cols: list[str] = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, lineterminator="\n", restval="")
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.9g}" if isinstance(v, float) else v) for k, v in r.items()})
    text = buf.getvalue()
    if path:
        Path(path).write_text(text, encoding="utf-8")
    return text
