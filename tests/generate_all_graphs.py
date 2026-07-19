#!/usr/bin/env python3
"""Generate the full dark-themed rockmoc graph suite into tests/graphs/.

Run once after install::

    py -3.12 -m pip install -e ".[plot]"
    py -3.12 tests/generate_all_graphs.py

All figures use the dark aerospace theme. Saving is optional in library APIs;
this script always writes PNGs under ``tests/graphs/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Prefer src/ layout when not installed
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib

matplotlib.use("Agg")  # headless; no GUI windows

from rockmoc.aerospike import design_spike
from rockmoc.analysis.compare import compare_nozzle_designs
from rockmoc.designers import DesignSpec, design_ideal, design_mln, design_tic
from rockmoc.designers.production import design_production
from rockmoc.designers.rao import design_rao
from rockmoc.plot import (
    area_ratio_sweep_rows,
    display_design_result,
    display_production_result,
    gamma_sweep_rows,
    of_ratio_sweep_rows,
    plot_aerospike,
    plot_area_ratio_sweep,
    plot_contour_overlay,
    plot_design_comparison,
    plot_gamma_sweep,
    plot_of_ratio_sweep,
    plot_performance_vs_altitude,
)
from rockmoc.plot.sweeps import altitude_sweep_rows
from rockmoc.thermo.perfect_gas import PerfectGas

OUT = Path(__file__).resolve().parent / "graphs"
OUT.mkdir(parents=True, exist_ok=True)

# Moderate mesh for reasonably fast generation
N_CHARS = 48
Rt = 0.05
Me = 3.0
gamma = 1.25
pc = 5.0e6
Tc = 3200.0


def main() -> list[str]:
    saved: list[str] = []
    gas = PerfectGas(gamma, R=350.0)
    spec = DesignSpec(
        gas=gas,
        throat_radius=Rt,
        exit_mach=Me,
        n_characteristics=N_CHARS,
    )

    print("== MLN design (with MOC net) ==")
    mln = design_mln(spec, store_net=True)
    figs = display_design_result(
        mln, show=False, save_dir=str(OUT), prefix="mln", show_moc=True
    )
    print(f"  MLN dashboard: {len(figs)} figures")
    for name in (
        "mln_contour.png",
        "mln_moc_field.png",
        "mln_wall_profiles.png",
        "mln_mach_pressure.png",
    ):
        p = OUT / name
        if p.exists():
            saved.append(str(p))

    print("== Ideal + TIC + Rao (contours) ==")
    ideal = design_ideal(spec, store_net=True)
    tic = design_tic(
        DesignSpec(
            gas=gas,
            throat_radius=Rt,
            exit_mach=Me,
            n_characteristics=N_CHARS,
            truncate_fraction=0.80,
        ),
        store_net=True,
    )
    rao = design_rao(spec, length_fraction=0.80, store_net=False)

    for res, prefix in ((ideal, "ideal"), (tic, "tic"), (rao, "rao")):
        display_design_result(
            res,
            show=False,
            save_dir=str(OUT),
            prefix=prefix,
            show_moc=res.net is not None,
        )
        for suffix in ("contour.png", "moc_field.png", "wall_profiles.png", "mach_pressure.png"):
            p = OUT / f"{prefix}_{suffix}"
            if p.exists():
                saved.append(str(p))

    print("== Contour overlay (family) ==")
    plot_contour_overlay(
        {
            "MLN": mln.contour,
            "Ideal": ideal.contour,
            "TIC 80%": tic.contour,
            "Rao TOP 80%": rao.contour,
        },
        title="Nozzle family overlay — MLN / Ideal / TIC / Rao",
        show=False,
        save_path=str(OUT / "family_contour_overlay.png"),
    )
    saved.append(str(OUT / "family_contour_overlay.png"))

    print("== Design comparison bars ==")
    cmp = compare_nozzle_designs(
        DesignSpec(gas=gas, throat_radius=Rt, exit_mach=Me, n_characteristics=N_CHARS),
        pc=pc,
        Tc=Tc,
    )
    plot_design_comparison(
        cmp.rows,
        title="Design family comparison (Me=3, γ=1.25)",
        show=False,
        save_path=str(OUT / "family_design_comparison.png"),
    )
    saved.append(str(OUT / "family_design_comparison.png"))

    print("== Performance vs altitude ==")
    alt_rows = altitude_sweep_rows(
        mln.contour,
        gas,
        throat_radius=Rt,
        exit_mach=Me,
        pc=pc,
        Tc=Tc,
    )
    plot_performance_vs_altitude(
        alt_rows,
        title="MLN performance vs altitude",
        show=False,
        save_path=str(OUT / "mln_performance_vs_altitude.png"),
    )
    saved.append(str(OUT / "mln_performance_vs_altitude.png"))

    print("== OF ratio sweep ==")
    of_rows = of_ratio_sweep_rows(
        throat_radius=Rt,
        exit_mach=Me,
        n_chars=24,
        pc=pc,
        Tc=Tc,
        use_cea=True,
    )
    plot_of_ratio_sweep(
        of_rows,
        title="O/F ratio sweep (LOX/RP-1 style)",
        show=False,
        save_path=str(OUT / "sweep_of_ratio.png"),
    )
    saved.append(str(OUT / "sweep_of_ratio.png"))

    print("== Gamma sweep ==")
    g_rows = gamma_sweep_rows(
        throat_radius=Rt,
        exit_mach=Me,
        n_chars=24,
        pc=pc,
        Tc=Tc,
    )
    plot_gamma_sweep(
        g_rows,
        title="Specific-heat ratio (γ) sweep",
        show=False,
        save_path=str(OUT / "sweep_gamma.png"),
    )
    saved.append(str(OUT / "sweep_gamma.png"))

    print("== Area-ratio sweep ==")
    eps_rows = area_ratio_sweep_rows(
        throat_radius=Rt,
        gamma=gamma,
        n_chars=24,
        pc=pc,
        Tc=Tc,
    )
    plot_area_ratio_sweep(
        eps_rows,
        title="Area-ratio (ε) design sweep",
        show=False,
        save_path=str(OUT / "sweep_area_ratio.png"),
    )
    saved.append(str(OUT / "sweep_area_ratio.png"))

    print("== Production pipeline dashboard ==")
    prod = design_production(
        DesignSpec(
            gas=gas,
            throat_radius=Rt,
            exit_mach=Me,
            n_characteristics=N_CHARS,
        ),
        kind="mln",
        apply_bl=True,
        pc=pc,
        Tc=Tc,
        store_net=True,
        n_chars_min=N_CHARS,  # keep demo generation fast
        certificate=False,
    )
    display_production_result(
        prod, show=False, save_dir=str(OUT), prefix="production"
    )
    for name in (
        "production_contour.png",
        "production_moc_field.png",
        "production_wall_profiles.png",
        "production_mach_pressure.png",
        "production_bl.png",
        "production_altitude.png",
    ):
        p = OUT / name
        if p.exists():
            saved.append(str(p))

    print("== Aerospike ==")
    spike = design_spike(
        Pc=pc,
        Tc=Tc,
        gamma=gamma,
        molar_mass=22.0,
        Re=0.15,
        expansion_ratio=25.0,
        Pa=101325.0,
        n_points=2000,
    )
    plot_aerospike(
        spike.X,
        spike.R,
        Re=spike.Re,
        title="Annular aerospike (Lee / Vernacchia)",
        show=False,
        save_path=str(OUT / "aerospike_contour.png"),
    )
    saved.append(str(OUT / "aerospike_contour.png"))

    # Deduplicate while preserving order
    seen: set[str] = set()
    unique: list[str] = []
    for p in saved:
        if p not in seen and Path(p).exists():
            seen.add(p)
            unique.append(p)

    manifest = OUT / "MANIFEST.txt"
    with open(manifest, "w", encoding="utf-8") as f:
        f.write("rockmoc dark-theme graph suite\n")
        f.write(f"output dir: {OUT}\n")
        f.write(f"count: {len(unique)}\n\n")
        for p in unique:
            f.write(Path(p).name + "\n")
    unique.append(str(manifest))

    print("\n========== GRAPHS WRITTEN ==========")
    for p in unique:
        print(f"  {p}")
    print(f"Total: {len(unique) - 1} PNGs + MANIFEST.txt")
    return unique


if __name__ == "__main__":
    main()
