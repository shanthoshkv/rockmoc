"""Generate student tutorial notebooks (run once if regenerating)."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def write_nb(name: str, cells: list[tuple[str, str]]) -> None:
    out_cells = []
    for typ, src in cells:
        raw = src.strip("\n")
        if not raw:
            lines: list[str] = []
        else:
            parts = raw.split("\n")
            lines = [p + "\n" for p in parts[:-1]] + [parts[-1]]
        if typ == "md":
            out_cells.append({"cell_type": "markdown", "metadata": {}, "source": lines})
        else:
            out_cells.append(
                {
                    "cell_type": "code",
                    "execution_count": None,
                    "metadata": {},
                    "outputs": [],
                    "source": lines,
                }
            )
    doc = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "cells": out_cells,
    }
    path = ROOT / name
    path.write_text(json.dumps(doc, indent=1), encoding="utf-8")
    print("wrote", path)


def main() -> None:
    write_nb(
        "01_first_mln.ipynb",
        [
            (
                "md",
                """# 01 — First MLN nozzle

**Goal:** Design an axisymmetric minimum-length nozzle, print a student report, export the wall, and (optionally) plot.

Install from the repo root:

```bash
pip install -e ".[plot]"
```
""",
            ),
            (
                "code",
                """from rockmoc import DesignSpec, design_mln, explainer, list_example_gases

print(explainer("mln"))
print()
print("Example gases:", list_example_gases())
""",
            ),
            (
                "code",
                """# Student-friendly spec: Rt in millimeters, modest mesh
spec = DesignSpec.student(Me=3.0, gamma=1.25, Rt_mm=50, n=40)
# Or: DesignSpec.student(Me=3.0, gas_name="lox_rp1", Rt_mm=50)

mln = design_mln(spec, store_net=True)
print(mln.report(length_unit="mm", angle_unit="deg"))
""",
            ),
            (
                "code",
                """# Export + lab report skeleton
mln.export_csv("nb01_mln_wall.csv")
path = mln.lab_report(
    "nb01_lab_report.md",
    student_name="Your Name",
    course="Propulsion Lab",
)
print("Wrote", path)
print("Summary keys:", list(mln.summary().keys())[:8])
""",
            ),
            (
                "code",
                """# Optional dark plots (needs matplotlib)
try:
    mln.plot(show=False, save_dir="nb01_figures", prefix="mln")
    print("Figures saved under nb01_figures/")
except ImportError as e:
    print("Plotting skipped:", e)
""",
            ),
            (
                "md",
                """## Checkpoint

- Is exit wall angle near 0 deg?
- Does area ratio look consistent with A/A*(Me) for gamma = 1.25?
- Read the **Limits** line in the report — what must CFD still do?
""",
            ),
        ],
    )

    write_nb(
        "02_compare_nozzles.ipynb",
        [
            (
                "md",
                """# 02 — Compare nozzle families

**Goal:** At fixed Me and gamma, compare MLN, ideal, TIC, Rao TOP, and a 15 deg cone.
Export a Markdown/CSV table for your report.
""",
            ),
            (
                "code",
                """from rockmoc import DesignSpec, compare_nozzle_designs, explainer

for k in ("mln", "ideal", "tic", "rao", "conical"):
    print(explainer(k))
    print()
""",
            ),
            (
                "code",
                """spec = DesignSpec.student(Me=3.0, gamma=1.25, Rt_mm=50, n=32)
cmp = compare_nozzle_designs(
    spec,
    rao_length_fraction=0.80,
    tic_fraction=0.80,
    cone_half_angle_deg=15.0,
    include_conical=True,
)
print(cmp.summary_table())
print()
print(cmp.to_markdown())
cmp.to_csv("nb02_compare.csv")
print("Wrote nb02_compare.csv")
""",
            ),
            (
                "md",
                """## Discussion

1. Which design is shortest in L/Rt?
2. Which has the largest exit angle? How does that hurt divergence efficiency?
3. Why is the cone still useful as a baseline?
""",
            ),
        ],
    )

    write_nb(
        "03_performance_altitude.ipynb",
        [
            (
                "md",
                """# 03 — Performance vs altitude

**Goal:** After designing an MLN, tabulate Cf and Isp from sea level to high altitude.
""",
            ),
            (
                "code",
                """from rockmoc import DesignSpec, design_mln, example_gas, off_design_altitude_sweep

gas = example_gas("lox_rp1")
spec = DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=40)
mln = design_mln(spec)
print(mln.report())

rows = off_design_altitude_sweep(
    mln.contour,
    gas,
    throat_radius=0.05,
    exit_mach=3.0,
    pc=5e6,
    Tc=3200.0,
    altitudes_m=[0, 5e3, 10e3, 20e3, 40e3, 80e3],
)
print("h_km, Cf, Isp_s")
for r in rows:
    print(r["altitude_m"] / 1000.0, r["Cf"], r["Isp_s"])
""",
            ),
            (
                "code",
                """try:
    from rockmoc.plot import plot_performance_vs_altitude

    plot_performance_vs_altitude(rows, show=False, save_path="nb03_altitude.png")
    print("Saved nb03_altitude.png")
except ImportError as e:
    print("Plot skipped:", e)
""",
            ),
            (
                "md",
                """## Discussion

Why does performance improve as ambient pressure falls?
When would a larger area ratio help vacuum but hurt sea-level thrust?
""",
            ),
        ],
    )

    write_nb(
        "04_nasa_planar.ipynb",
        [
            (
                "md",
                """# 04 — NASA planar MLN (TM-X-1502)

**Goal:** Reproduce planar design rules: theta_max = nu_e/2 and y_e/y_t = A/A*.
""",
            ),
            (
                "code",
                """import math
from rockmoc import PerfectGas, design_planar_mln_nasa, explainer

print(explainer("planar_mln"))

# Paper sample-style case (coarser dnu for notebook speed)
Me, g = 4.05, 1.36
gas = PerfectGas(g)
nu_e = float(gas.prandtl_meyer(Me))
r = design_planar_mln_nasa(Me, g, dnu_deg=0.5, throat_half_height=1.0)
print("theta_max_deg", math.degrees(r.theta_max_rad))
print("nu_e/2_deg", math.degrees(0.5 * nu_e))
print("y_e", r.contour.r[-1], "A/A*", float(gas.area_ratio(Me)))
print("x_e_unit_ht", r.contour.x[-1])
""",
            ),
            (
                "code",
                """# Classic gamma=1.4, Me=3 teaching case
r2 = design_planar_mln_nasa(3.0, 1.4, dnu_deg=1.0, throat_half_height=0.05)
print(r2.summary())
r2.contour.export_csv("nb04_planar_wall.csv")
""",
            ),
            (
                "md",
                """## Checkpoint

Match theta_max to nu_e/2. For unit throat half-height, y_e must equal A/A*(Me).
""",
            ),
        ],
    )


if __name__ == "__main__":
    main()
