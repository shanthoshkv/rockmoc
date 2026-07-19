"""Command-line interface for rockmoc."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

from rockmoc._version import __version__
from rockmoc.designers import DesignSpec, design_ideal, design_mln, design_tic
from rockmoc.thermo import PerfectGas

try:
    import typer
except ImportError:  # pragma: no cover
    typer = None  # type: ignore


def _run_design(
    kind: str,
    gamma: float,
    R: float,
    Rt: float,
    Me: Optional[float],
    area_ratio: Optional[float],
    n: int,
    truncate: float,
    output: Path,
    json_meta: Optional[Path],
) -> None:
    gas = PerfectGas(gamma=gamma, R=R)
    spec = DesignSpec(
        gas=gas,
        throat_radius=Rt,
        exit_mach=Me,
        area_ratio=area_ratio,
        n_characteristics=n,
        truncate_fraction=truncate if kind == "tic" else None,
    )
    designers = {"mln": design_mln, "ideal": design_ideal, "tic": design_tic}
    if kind not in designers:
        raise SystemExit(f"Unknown design kind {kind!r}; choose mln|ideal|tic")
    result = designers[kind](spec)
    result.export_csv(str(output))
    summary = result.summary()
    print(json.dumps(summary, indent=2))
    if json_meta is not None:
        json_meta.write_text(json.dumps(summary, indent=2), encoding="utf-8")


if typer is not None:
    app = typer.Typer(
        name="rockmoc",
        help="Axisymmetric MOC rocket nozzle design tool",
        add_completion=False,
    )

    @app.callback()
    def main(
        version: bool = typer.Option(
            False, "--version", help="Show version and exit."
        ),
    ) -> None:
        if version:
            typer.echo(f"rockmoc {__version__}")
            raise typer.Exit()

    @app.command("design")
    def design_cmd(
        type: str = typer.Option("mln", "--type", "-t", help="mln | ideal | tic"),
        gamma: float = typer.Option(1.4, "--gamma", "-g"),
        R: float = typer.Option(287.05, "--R", help="Specific gas constant [J/(kg·K)]"),
        Rt: float = typer.Option(0.05, "--Rt", help="Throat radius [m]"),
        Me: Optional[float] = typer.Option(None, "--Me", help="Design exit Mach"),
        area_ratio: Optional[float] = typer.Option(
            None, "--area-ratio", "-e", help="Design A_e/A*"
        ),
        n: int = typer.Option(40, "--n", help="Number of characteristics"),
        truncate: float = typer.Option(
            0.8, "--truncate", help="TIC length fraction of ideal"
        ),
        output: Path = typer.Option(Path("nozzle_wall.csv"), "--output", "-o"),
        json_meta: Optional[Path] = typer.Option(
            None, "--json", help="Optional path for JSON summary"
        ),
    ) -> None:
        """Design a nozzle contour and write wall CSV."""
        if Me is None and area_ratio is None:
            Me = 2.5
        _run_design(
            type.lower(), gamma, R, Rt, Me, area_ratio, n, truncate, output, json_meta
        )

else:
    # Minimal fallback without typer
    def app(argv: list[str] | None = None) -> None:  # type: ignore[misc]
        argv = list(sys.argv[1:] if argv is None else argv)
        if not argv or argv[0] in {"-h", "--help"}:
            print(
                "Usage: rockmoc design --type mln|ideal|tic --Me 2.5 --Rt 0.05 -o wall.csv\n"
                "Install typer for the full CLI: pip install rockmoc[cli]"
            )
            return
        if argv[0] == "--version":
            print(f"rockmoc {__version__}")
            return
        # very small parser
        kind = "mln"
        gamma = 1.4
        R = 287.05
        Rt = 0.05
        Me: float | None = 2.5
        area_ratio: float | None = None
        n = 40
        truncate = 0.8
        output = Path("nozzle_wall.csv")
        json_meta = None
        args = argv[1:] if argv[0] == "design" else argv
        i = 0
        while i < len(args):
            a = args[i]
            if a in {"--type", "-t"}:
                kind = args[i + 1]
                i += 2
            elif a in {"--gamma", "-g"}:
                gamma = float(args[i + 1])
                i += 2
            elif a == "--R":
                R = float(args[i + 1])
                i += 2
            elif a == "--Rt":
                Rt = float(args[i + 1])
                i += 2
            elif a == "--Me":
                Me = float(args[i + 1])
                area_ratio = None
                i += 2
            elif a in {"--area-ratio", "-e"}:
                area_ratio = float(args[i + 1])
                Me = None
                i += 2
            elif a == "--n":
                n = int(args[i + 1])
                i += 2
            elif a == "--truncate":
                truncate = float(args[i + 1])
                i += 2
            elif a in {"--output", "-o"}:
                output = Path(args[i + 1])
                i += 2
            elif a == "--json":
                json_meta = Path(args[i + 1])
                i += 2
            else:
                print(f"Unknown arg: {a}", file=sys.stderr)
                raise SystemExit(2)
        _run_design(kind, gamma, R, Rt, Me, area_ratio, n, truncate, output, json_meta)


def main() -> None:
    if typer is not None:
        app()
    else:
        app()  # type: ignore[call-arg]


if __name__ == "__main__":
    main()
