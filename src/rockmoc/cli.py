"""Command-line interface (standard-library ``argparse``; no extra dependency).

Examples::

    rockmoc design --type ideal --Me 3 --gamma 1.2 --Rt 0.05 --export csv,dxf --units mm
    rockmoc design --Me 3 --Rt 50 --in-units mm --units in --points 200 --export csv,txt
    rockmoc design --area-ratio 25 --propellants LOX/RP1 --of 2.4 --pc 7e6  # real gas + thrust
    rockmoc sweep --propellants LOX/RP1 --of 2.4 --pc 7e6 --vary gas_model=frozen,equilibrium
    rockmoc design --Me 3 --stdout > wall.csv                                # wall CSV to stdout
    rockmoc sweep --vary Me=2:4:5 --vary of=2.0:3.0:6 --propellants LOX/RP1 --pc 7e6 --out s.csv
    rockmoc sweep --vary pc_bar=30:100:8 --vary pe_pc=0.005,0.01 --jobs 4 --export csv -o cases
    rockmoc run nozzle.toml          # a [sweep] table in the TOML runs a sweep
    rockmoc validate nozzle.toml
    rockmoc compare a.csv b.csv --units mm
    rockmoc --version
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from rockmoc._version import __version__
from rockmoc.export import UNITS
from rockmoc.spec_io import (
    FORMATS,
    KINDS,
    build_design,
    export_design,
    load_spec,
    performance_from_spec,
    validate_spec,
)


def _base() -> argparse.ArgumentParser:
    """Design inputs shared by ``design`` and ``sweep``."""
    b = argparse.ArgumentParser(add_help=False)
    b.add_argument("-t", "--type", choices=KINDS, default="mln")
    b.add_argument("-g", "--gamma", type=float, default=1.4)
    b.add_argument("--R", type=float, default=287.05, help="gas constant [J/(kg K)]")
    b.add_argument("--propellants", metavar="OX/FUEL",
                   help="e.g. LOX/RP1: real-gas RocketCEA table at --of, --pc (overrides -g/--R)")
    b.add_argument("--gas-model", choices=("perfect", "equilibrium", "frozen"),
                   help="perfect = constant γ; equilibrium|frozen = real gas, γ(T) varies along "
                        "the nozzle (default: equilibrium with --propellants, else perfect)")
    b.add_argument("--eps-max", type=float, help="real-gas table limit in area ratio (default 400)")
    b.add_argument("--of", type=float, help="oxidizer/fuel mass ratio (with --propellants)")
    b.add_argument("--Rt", type=float, default=0.05, help="throat radius [--in-units]")
    b.add_argument("--in-units", choices=UNITS, default="m", help="units of --Rt (default m)")
    exp = b.add_mutually_exclusive_group()
    exp.add_argument("--Me", type=float, help="design exit Mach")
    exp.add_argument("-e", "--area-ratio", type=float, help="design A_e/A*")
    exp.add_argument("--pe-pc", type=float, help="design exit/chamber pressure ratio")
    b.add_argument("-n", "--n", type=int, default=80, help="number of characteristics (CFD: axis Mach rms <1 % needs n >= 80)")
    b.add_argument("--rc", type=float, default=1.5, help="throat Rc/Rt (ideal/tic)")
    b.add_argument("--truncate", type=float, default=0.8, help="TIC length fraction of ideal")
    b.add_argument("--length", type=float, default=0.8, help="Rao length fraction of 15° cone")
    b.add_argument("--pc", type=float, help="chamber pressure [Pa] → adds thrust/Isp/mdot")
    b.add_argument("--Tc", type=float, help="chamber temperature [K] (default: CEA, else 3200)")
    b.add_argument("--pa", type=float, help="ambient pressure [Pa]")
    b.add_argument("--altitude", type=float, help="ambient from altitude [m] (std atmosphere)")
    return b


def _formats(s: str) -> str:
    bad = [f for f in (x.strip() for x in s.split(",")) if f and f not in FORMATS]
    if bad:
        raise argparse.ArgumentTypeError(f"unknown formats {bad}; choose from {FORMATS}")
    return s


def _out_args(p: argparse.ArgumentParser, default_export: str) -> None:
    p.add_argument("--export", type=_formats, default=default_export, help=f"comma list of {','.join(FORMATS)}")
    p.add_argument("--units", choices=UNITS, default="mm", help="output units (default mm)")
    p.add_argument("--points", type=int, help="resample exported wall to N points")
    p.add_argument("--thickness", type=float, help="outer-wall offset for 'shell' [--units]")
    p.add_argument("-o", "--out-dir", default=".")
    p.add_argument("--name", default="nozzle")


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="rockmoc", description="MOC rocket nozzle design.")
    ap.add_argument("-V", "--version", action="version", version=f"rockmoc {__version__}")
    sub = ap.add_subparsers(dest="cmd")

    d = sub.add_parser("design", parents=[_base()], help="design a nozzle from CLI options")
    _out_args(d, "csv")
    d.add_argument("--stdout", action="store_true", help="print the wall CSV to stdout")

    r = sub.add_parser("run", help="design (or sweep, if it has [sweep]) from a TOML spec")
    r.add_argument("spec", type=Path)
    r.add_argument("--stdout", action="store_true", help="print the wall CSV to stdout")

    s = sub.add_parser("sweep", parents=[_base()], help="Cartesian parameter sweep → summary CSV")
    s.add_argument("spec", nargs="?", type=Path, help="optional TOML base spec (else flags)")
    s.add_argument("--vary", action="append", default=[], metavar="NAME=a:b:n|a,b,c",
                   help="repeatable. NAME: Me area_ratio pe_pc pe pc Tc pa altitude of gamma R "
                        "Rt n rc type length truncate (pressures accept _bar/_mpa/_kpa/_psia)")
    s.add_argument("--out", type=Path, help="summary CSV path (default: stdout)")
    s.add_argument("--jobs", type=int, default=1, help="parallel worker processes")
    _out_args(s, "")

    v = sub.add_parser("validate", help="check a TOML spec without running the solver")
    v.add_argument("spec", type=Path)

    c = sub.add_parser("compare", help="deviation between two wall CSV/TXT files")
    c.add_argument("a", type=Path)
    c.add_argument("b", type=Path)
    c.add_argument("--units", choices=UNITS, default="mm", help="units of both files")
    return ap


def _spec_from_args(a: argparse.Namespace) -> dict:
    design = {"type": a.type, "throat_radius": a.Rt, "units": a.in_units,
              "n_characteristics": a.n, "throat_rc_over_rt": a.rc,
              "truncate_fraction": a.truncate, "length_fraction": a.length}
    if a.area_ratio is not None:
        design["area_ratio"] = a.area_ratio
    elif a.pe_pc is not None:
        design["pressure_ratio"] = a.pe_pc
    else:
        design["exit_mach"] = a.Me if a.Me is not None else 2.5
    spec: dict = {"gas": {"gamma": a.gamma, "R": a.R}, "design": design,
                  "export": {"formats": [f.strip() for f in a.export.split(",") if f.strip()],
                             "units": a.units, "points": a.points, "out_dir": a.out_dir,
                             "name": a.name}}
    if a.thickness is not None:
        spec["export"]["thickness"] = a.thickness
    perf = {k: v for k, v in (("pc", a.pc), ("Tc", a.Tc), ("pa", a.pa),
                               ("altitude", a.altitude)) if v is not None}
    if perf:
        spec["performance"] = perf
    if a.gas_model:
        spec["gas"]["model"] = a.gas_model
    if a.eps_max:
        spec["gas"]["eps_max"] = a.eps_max
    if a.propellants:
        if a.of is None:
            raise SystemExit("--propellants needs --of")
        ox, _, fuel = a.propellants.partition("/")
        spec["propellants"] = {"ox": ox, "fuel": fuel or "RP1", "of": a.of}
    return spec


def _axes_from_toml(table: dict) -> dict:
    from rockmoc.sweep import parse_axis

    return {k: (parse_axis(f"{k}={v}")[1] if isinstance(v, str) else list(v))
            for k, v in table.items()}


def _sweep(spec: dict, axes: dict, a: argparse.Namespace | None, out: Path | None) -> None:
    from rockmoc.sweep import run_sweep, write_rows

    fmts = [f.strip() for f in (a.export if a else "").split(",") if f.strip()]
    rows = run_sweep(spec, axes, export_dir=a.out_dir if fmts else None, export_formats=fmts,
                     jobs=a.jobs if a else 1)
    text = write_rows(rows, out)
    if out is None:
        sys.stdout.write(text)
    bad = [r for r in rows if "error" in r]
    print(f"{len(rows) - len(bad)}/{len(rows)} cases ok"
          + (f" → {out}" if out else ""), file=sys.stderr)
    for r in bad[:5]:
        print(f"  case {r['case']}: {r['error']}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    """Console-script entry point."""
    args = _parser().parse_args(argv)
    if args.cmd is None:
        _parser().print_help()
        return 0
    if args.cmd == "compare":
        from rockmoc.analysis.compare import compare_contours
        from rockmoc.geometry.contour import Contour

        a, b = (Contour.from_csv(p, units=args.units) for p in (args.a, args.b))
        print(json.dumps(compare_contours(a, b), indent=2))
        return 0
    if args.cmd == "sweep":
        from rockmoc.sweep import parse_axis

        spec = load_spec(args.spec) if args.spec else _spec_from_args(args)
        axes = dict(parse_axis(v) for v in args.vary)
        axes.update(_axes_from_toml(spec.pop("sweep", {})))
        _sweep(spec, axes, args, args.out)
        return 0
    spec = _spec_from_args(args) if args.cmd == "design" else load_spec(args.spec)
    if args.cmd == "validate":
        validate_spec(spec)
        print(f"{args.spec}: OK")
        return 0
    if args.cmd == "run" and "sweep" in spec:
        axes = _axes_from_toml(spec.pop("sweep"))
        _sweep(spec, axes, None, Path(spec.get("export", {}).get("sweep_csv", "sweep.csv")))
        return 0
    validate_spec(spec)
    result = build_design(spec)
    paths = export_design(result, spec)
    summary = result.summary()
    perf = performance_from_spec(result, spec)
    if perf:
        summary["performance"] = perf
    if args.stdout:
        import numpy as np

        c = result.contour
        np.savetxt(sys.stdout, np.column_stack([c.x, c.r]), delimiter=",",
                   header="x_m,r_m", comments="", fmt="%.9g")
        print(json.dumps(summary, indent=2, default=str), file=sys.stderr)
    else:
        print(json.dumps(summary, indent=2, default=str))
    for p in paths:
        print(f"wrote {p}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
