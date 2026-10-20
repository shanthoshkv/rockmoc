"""TOML design specifications — reproducible, version-controllable nozzle designs.

Example ``nozzle.toml``::

    [gas]
    gamma = 1.2
    R = 350.0                 # or: molecular_weight = 23.7

    [design]
    type = "ideal"            # mln | ideal | tic | rao
    throat_radius = 0.05      # m
    exit_mach = 3.0           # or: area_ratio = 4.9
    n_characteristics = 60
    throat_rc_over_rt = 1.5
    truncate_fraction = 0.8   # tic
    length_fraction = 0.8     # rao (fraction of 15° cone)

    [boundary_layer]          # optional δ* wall offset
    pc = 5.0e6
    Tc = 3300.0

    [export]
    formats = ["csv", "dxf", "openfoam"]
    units = "mm"
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
FORMATS = ("csv", "cad", "dxf", "openfoam", "json")


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


def build_design(spec: dict[str, Any]) -> Any:
    """Run the design described by a spec dict; returns a DesignResult."""
    from rockmoc.designers import DesignSpec, design_ideal, design_mln, design_tic
    from rockmoc.designers.rao_control import design_rao_control
    from rockmoc.thermo import PerfectGas

    g = dict(spec.get("gas", {}))
    if "molecular_weight" in g:
        gas = PerfectGas.from_molecular_weight(g.get("gamma", 1.4), g["molecular_weight"])
    else:
        gas = PerfectGas(gamma=g.get("gamma", 1.4), R=g.get("R", 287.052874))
    d = dict(spec.get("design", {}))
    kind = d.pop("type", "mln")
    length_fraction = d.pop("length_fraction", 0.8)
    if kind != "tic":
        d.pop("truncate_fraction", None)
    elif "truncate_fraction" not in d:
        d["truncate_fraction"] = 0.8
    ds = DesignSpec(gas=gas, **d)
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
    units = e.get("units", "mm")
    written: list[Path] = []
    for fmt in e.get("formats", ["csv"]):
        if fmt == "csv":
            written.append(ex.export_csv(result.contour, out / f"{name}.csv", units=units))
        elif fmt == "cad":
            written.append(ex.export_cad_points(result.contour, out / f"{name}_xyz.txt", units=units))
        elif fmt == "dxf":
            written.append(ex.export_dxf(result.contour, out / f"{name}.dxf", units=units,
                                         mirror=e.get("mirror", False)))
        elif fmt == "openfoam":
            of = e.get("openfoam", {})
            written.append(ex.export_blockmesh(result.contour, out / name / "system" / "blockMeshDict", **of))
        elif fmt == "json":
            p = out / f"{name}.json"
            p.write_text(json.dumps(result.summary(), indent=2, default=str), encoding="utf-8")
            written.append(p)
    return written
