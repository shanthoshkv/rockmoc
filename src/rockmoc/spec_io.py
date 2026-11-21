"""TOML design specifications — reproducible, version-controllable nozzle designs.

Example ``nozzle.toml``::

    [gas]
    gamma = 1.2
    R = 350.0                 # or: molecular_weight = 23.7
    # model = "perfect"       # perfect (constant γ; default without [propellants]) |
    #                         # equilibrium (default with [propellants]) | frozen |
    #                         # multi_gamma (gamma_c, gamma_t, gamma_e, R)
    # eps_max = 400           # real-gas table limit (area ratio); larger designs are rejected

    [design]
    type = "ideal"            # mln | ideal | tic | rao
    throat_radius = 0.05      # in [design] units (default m)
    units = "m"               # input length units: m | cm | mm | in | ft
    exit_mach = 3.0           # or: area_ratio = 4.9
    n_characteristics = 80
    throat_rc_over_rt = 1.5
    truncate_fraction = 0.8   # tic
    length_fraction = 0.8     # rao (fraction of 15° cone)

    [propellants]             # RocketCEA real-gas table: γ(T), a, ρ, ν(M) vary along the nozzle
    ox = "LOX"
    fuel = "RP1"
    of = 2.4                  # oxidizer/fuel mass ratio
    pc = 5.0e6                # chamber pressure [Pa] (default: [performance] pc)

    [performance]             # optional: thrust / Isp / mdot in the printed summary
    pc = 5.0e6                # chamber pressure [Pa]
    Tc = 3300.0               # K (default: from [propellants], else 3200)
    pa = 0.0                  # ambient [Pa]   (or: altitude = 10000  # m)

    [boundary_layer]          # optional δ* wall offset
    pc = 5.0e6
    Tc = 3300.0

    [export]
    formats = ["csv", "dxf", "openfoam"]
    units = "mm"              # output units (default: [design] units)
    points = 200              # optional: resample wall to N points
    thickness = 3.0           # "shell" format: outer wall offset, in export units
    out_dir = "out"
    name = "nozzle"
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

try:
    import tomllib
except ImportError:  # Python 3.10
    import tomli as tomllib  # type: ignore[no-redef]

from rockmoc.exceptions import DesignError

KINDS = ("mln", "ideal", "tic", "rao")
FORMATS = ("csv", "txt", "cad", "sldcrv", "shell", "dxf", "openfoam", "json")


def load_spec(path: str | Path) -> dict[str, Any]:
    """Read and validate a TOML spec; returns the parsed dictionary."""
    with open(path, "rb") as fh:
        spec = tomllib.load(fh)
    design = spec.get("design", {})
    if design.get("type", "mln") not in KINDS:
        raise DesignError(f"[design] type must be one of {KINDS}")
    bad = [f for f in spec.get("export", {}).get("formats", []) if f not in FORMATS]
    if bad:
        raise DesignError(f"[export] unknown formats {bad}; choose from {FORMATS}")
    return spec


GAS_MODELS = ("perfect", "equilibrium", "frozen", "multi_gamma")


def resolve_gas(spec: dict[str, Any]) -> tuple[Any, float | None]:
    """``(gas, Tc)`` for a spec. ``Tc`` (K) may be None.

    ``[gas] model``: ``perfect`` (constant γ), ``equilibrium`` / ``frozen`` (real-gas
    table from RocketCEA at ``[propellants]`` — γ, T, ρ, a, ν(M) all vary along the
    nozzle), ``multi_gamma`` (γ_c/γ_t/γ_e blended). Default is ``equilibrium`` when
    ``[propellants]`` is given (falls back to ``perfect`` with a warning if RocketCEA is
    missing), else ``perfect``.
    """
    import warnings

    from rockmoc.thermo import PerfectGas, VariableGammaGas

    g = dict(spec.get("gas", {}))
    pr = spec.get("propellants")
    perf = spec.get("performance", {})
    explicit = "model" in g
    model = g.get("model", "equilibrium" if pr else "perfect")
    if model not in GAS_MODELS:
        raise DesignError(f"[gas] model must be one of {GAS_MODELS}; got {model!r}")
    if model in ("equilibrium", "frozen"):
        if not pr:
            raise DesignError(f"[gas] model = {model!r} needs [propellants] (ox, fuel, of)")
        pc = pr.get("pc", perf.get("pc", 5e6))
        try:
            from rockmoc.thermo.cea_table import build_cea_table

            eps_req = spec.get("design", {}).get("area_ratio", 0)
            eps_max = g.get("eps_max", max(400.0, 1.5 * eps_req))
            tb = build_cea_table(pr.get("ox", "LOX"), pr.get("fuel", "RP1"), pr["of"], pc,
                                 mode=model, eps_max=eps_max, n=g.get("table_points", 70))
            return VariableGammaGas(tb), tb.Tc
        except ImportError:
            if explicit:
                raise
            warnings.warn("RocketCEA not installed: using a constant-γ stand-in gas "
                          "(pip install rocketcea for real-gas tables)", stacklevel=2)
            model = "perfect"
    if model == "multi_gamma":
        from rockmoc.thermo import build_multi_gamma_table

        R = g.get("R", 287.052874)
        if "molecular_weight" in g:
            R = 8314.462618 / g["molecular_weight"]
        Tc = perf.get("Tc", 3200.0)
        tb = build_multi_gamma_table(gamma_c=g["gamma_c"], gamma_t=g["gamma_t"],
                                     gamma_e=g["gamma_e"], R=R, pc=perf.get("pc", 5e6), Tc=Tc,
                                     M_max=g.get("M_max", 8.0), n=g.get("table_points", 400))
        return VariableGammaGas(tb), Tc
    if pr:
        from rockmoc.thermo.propellant import of_props

        pc = pr.get("pc", perf.get("pc", 5e6))
        p = of_props(pr.get("ox", "LOX"), pr.get("fuel", "RP1"), pr["of"], pc)
        return PerfectGas(p["gamma"], p["R"]), p["Tc"]
    if "molecular_weight" in g:
        return PerfectGas.from_molecular_weight(g.get("gamma", 1.4), g["molecular_weight"]), None
    return PerfectGas(gamma=g.get("gamma", 1.4), R=g.get("R", 287.052874)), None


def _design_spec(spec: dict[str, Any]) -> tuple[Any, str, float]:
    from rockmoc.designers import DesignSpec

    gas, _ = resolve_gas(spec)
    d = dict(spec.get("design", {}))
    kind = d.pop("type", "mln")
    length_fraction = d.pop("length_fraction", 0.8)
    if kind != "tic":
        d.pop("truncate_fraction", None)
    elif "truncate_fraction" not in d:
        d["truncate_fraction"] = 0.8
    ds = DesignSpec.from_units(gas=gas, units=d.pop("units", "m"), **d)
    return ds, kind, length_fraction


def validate_spec(spec: dict[str, Any]) -> None:
    """Check a spec without running the solver; raises DesignError/ValueError."""
    from rockmoc.export import UNITS

    ds, _, _ = _design_spec(spec)
    ds.validate()
    e = spec.get("export", {})
    for u in (e.get("units"), spec.get("design", {}).get("units")):
        if u is not None and u not in UNITS:
            raise DesignError(f"units {u!r} not in {UNITS}")
    if e.get("points") is not None and int(e["points"]) < 2:
        raise DesignError("[export] points must be ≥ 2")
    p = spec.get("performance", {})
    if p.get("pc", 1.0) <= 0 or p.get("Tc", 1.0) <= 0 or p.get("pa", 0.0) < 0:
        raise DesignError("[performance] pc, Tc must be > 0 and pa ≥ 0")


def performance_from_spec(result: Any, spec: dict[str, Any]) -> dict | None:
    """Thrust/Isp/mdot summary if ``[performance]`` or ``[propellants]`` is present."""
    if "performance" not in spec and "propellants" not in spec:
        return None
    from rockmoc.analysis.atmosphere import pressure_at_altitude

    p = spec.get("performance", {})
    _, Tc_cea = resolve_gas(spec)
    pa = pressure_at_altitude(p["altitude"]) if "altitude" in p else p.get("pa", 0.0)
    return result.performance(pc=p.get("pc", spec.get("propellants", {}).get("pc", 5e6)),
                              Tc=p.get("Tc", Tc_cea or 3200.0), pa=pa)


def build_design(spec: dict[str, Any]) -> Any:
    """Run the design described by a spec dict; returns a DesignResult."""
    from rockmoc.designers import design_ideal, design_mln, design_tic
    from rockmoc.designers.rao_control import design_rao_control

    ds, kind, length_fraction = _design_spec(spec)
    if kind == "rao":
        result = design_rao_control(ds, length_fraction=length_fraction)
    else:
        result = {"mln": design_mln, "ideal": design_ideal, "tic": design_tic}[kind](ds)
    bl = spec.get("boundary_layer")
    if bl:
        result = result.with_boundary_layer(**bl)
    return result


def export_design(result: Any, spec: dict[str, Any]) -> list[Path]:
    """Write the formats listed in ``[export]``; returns the written paths."""
    from rockmoc import export as ex

    e = spec.get("export", {})
    out = Path(e.get("out_dir", "."))
    out.mkdir(parents=True, exist_ok=True)
    name = e.get("name", "nozzle")
    units = e.get("units", spec.get("design", {}).get("units", "mm"))
    n = e.get("points")
    written: list[Path] = []
    for fmt in e.get("formats", ["csv"]):
        if fmt == "csv":
            written.append(ex.export_csv(result.contour, out / f"{name}.csv", units=units, n=n))
        elif fmt == "txt":
            written.append(ex.export_contour(result.contour, out / f"{name}.txt", units=units, n=n))
        elif fmt == "cad":
            written.append(ex.export_cad_points(result.contour, out / f"{name}_xyz.txt",
                                                units=units, n=n))
        elif fmt == "sldcrv":  # SolidWorks Insert > Curve > Curve Through XYZ Points
            written.append(ex.export_cad_points(result.contour, out / f"{name}.sldcrv",
                                                units=units, n=n))
        elif fmt == "shell":
            written.append(ex.export_shell(result.contour, out / f"{name}_shell.csv",
                                           e.get("thickness", 3.0), units=units, n=n))
        elif fmt == "dxf":
            written.append(ex.export_dxf(result.contour, out / f"{name}.dxf", units=units,
                                         mirror=e.get("mirror", False)))
        elif fmt == "openfoam":
            of = e.get("openfoam", {})
            bm = out / name / "system" / "blockMeshDict"
            written.append(ex.export_blockmesh(result.contour, bm, **{"units": units, **of}))
        elif fmt == "json":
            p = out / f"{name}.json"
            doc = {"units": {"summary": "SI (m, Pa, K, rad)", "wall_export": units},
                   **result.summary()}
            perf = performance_from_spec(result, spec)
            if perf:
                doc["performance"] = perf
            p.write_text(json.dumps(doc, indent=2, default=str), encoding="utf-8")
            written.append(p)
    return written
