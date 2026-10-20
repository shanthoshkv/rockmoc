"""Command-line interface (standard-library ``argparse``; no extra dependency).

Examples::

    rockmoc design --type ideal --Me 3 --gamma 1.2 --Rt 0.05 --export csv,dxf --units mm
    rockmoc design --type rao --area-ratio 25 --length 0.8 --export openfoam
    rockmoc run nozzle.toml
    rockmoc --version
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from rockmoc._version import __version__
from rockmoc.spec_io import FORMATS, KINDS, build_design, export_design, load_spec


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="rockmoc", description="MOC rocket nozzle design.")
    ap.add_argument("-V", "--version", action="version", version=f"rockmoc {__version__}")
    sub = ap.add_subparsers(dest="cmd")

    d = sub.add_parser("design", help="design a nozzle from command-line options")
    d.add_argument("-t", "--type", choices=KINDS, default="mln")
    d.add_argument("-g", "--gamma", type=float, default=1.4)
    d.add_argument("--R", type=float, default=287.05, help="gas constant [J/(kg K)]")
    d.add_argument("--Rt", type=float, default=0.05, help="throat radius [m]")
    exp = d.add_mutually_exclusive_group()
    exp.add_argument("--Me", type=float, help="design exit Mach")
    exp.add_argument("-e", "--area-ratio", type=float, help="design A_e/A*")
    d.add_argument("-n", "--n", type=int, default=40, help="number of characteristics")
    d.add_argument("--rc", type=float, default=1.5, help="throat Rc/Rt (ideal/tic)")
    d.add_argument("--truncate", type=float, default=0.8, help="TIC length fraction of ideal")
    d.add_argument("--length", type=float, default=0.8, help="Rao length fraction of 15° cone")
    d.add_argument("--export", default="csv", help=f"comma list of {','.join(FORMATS)}")
    d.add_argument("--units", choices=("m", "mm", "in"), default="mm")
    d.add_argument("-o", "--out-dir", default=".")
    d.add_argument("--name", default="nozzle")

    r = sub.add_parser("run", help="design from a TOML spec file")
    r.add_argument("spec", type=Path)
    return ap


def _spec_from_args(a: argparse.Namespace) -> dict:
    design = {"type": a.type, "throat_radius": a.Rt, "n_characteristics": a.n,
              "throat_rc_over_rt": a.rc, "truncate_fraction": a.truncate,
              "length_fraction": a.length}
    if a.area_ratio is not None:
        design["area_ratio"] = a.area_ratio
    else:
        design["exit_mach"] = a.Me if a.Me is not None else 2.5
    return {"gas": {"gamma": a.gamma, "R": a.R}, "design": design,
            "export": {"formats": [f.strip() for f in a.export.split(",") if f.strip()],
                       "units": a.units, "out_dir": a.out_dir, "name": a.name}}


def main(argv: list[str] | None = None) -> int:
    """Console-script entry point."""
    args = _parser().parse_args(argv)
    if args.cmd is None:
        _parser().print_help()
        return 0
    spec = load_spec(args.spec) if args.cmd == "run" else _spec_from_args(args)
    result = build_design(spec)
    paths = export_design(result, spec)
    print(json.dumps(result.summary(), indent=2, default=str))
    for p in paths:
        print(f"wrote {p}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
