# How to use rockmoc

A friendly, detailed guide to install, design nozzles, plot results, export walls, and validate your setup.

If you only want a one-page overview, start with the [README](../README.md).  
If you want theory, see [theory.md](theory.md) and the NASA/aerospike notes.  
If you want validation methodology, see [validation_deep.md](validation_deep.md).

---

## Contents

1. [What you need](#1-what-you-need)
2. [Install (5 minutes)](#2-install-5-minutes)
3. [Mental model: how a design works](#3-mental-model-how-a-design-works)
4. [Your first nozzle (step by step)](#4-your-first-nozzle-step-by-step)
5. [Choosing inputs: Me vs area ratio vs pe/pc](#5-choosing-inputs-me-vs-area-ratio-vs-pepc)
6. [Which nozzle type should I use?](#6-which-nozzle-type-should-i-use)
7. [Recipe book (copy-paste)](#7-recipe-book-copy-paste)
8. [Reading the results](#8-reading-the-results)
9. [Student helpers (reports, gases, tables, cone)](#9-student-helpers-reports-gases-tables-cone)
10. [Plots and figures](#10-plots-and-figures)
11. [Export for CAD / CFD](#11-export-for-cad--cfd)
12. [Performance, altitude, and sweeps](#12-performance-altitude-and-sweeps)
13. [Production pipeline (recommended for “final” contours)](#13-production-pipeline-recommended-for-final-contours)
14. [Thermochemistry (CEA / Cantera / CoolProp)](#14-thermochemistry-cea--cantera--coolprop)
15. [Command-line usage](#15-command-line-usage)
16. [Mesh density: speed vs accuracy](#16-mesh-density-speed-vs-accuracy)
17. [Validation: prove the install](#17-validation-prove-the-install)
18. [Common mistakes & FAQ](#18-common-mistakes--faq)
19. [Where to go next](#19-where-to-go-next)

**Also for students:** [glossary.md](glossary.md) · [student_problems.md](student_problems.md) · [nozzle_explainers.md](nozzle_explainers.md) · [notebooks](../examples/notebooks/)

---

## 1. What you need

| Item | Notes |
|------|--------|
| **Python** | 3.10 or newer (3.12 is well tested) |
| **OS** | Windows, macOS, or Linux |
| **Core packages** | NumPy, SciPy (installed automatically) |
| **Optional plots** | Matplotlib via `rockmoc[plot]` |
| **CLI** | built in (`argparse`): `rockmoc design …`, `rockmoc run spec.toml` |
| **Optional chemistry** | RocketCEA, Cantera, or CoolProp |

You do **not** need CFD, MATLAB, or a GPU.

---

## 2. Install (5 minutes)

### From a clone (developers)

```bash
cd path/to/rockmoc

# Create a virtual environment (recommended)
py -3.12 -m venv .venv

# Windows PowerShell
.\.venv\Scripts\Activate.ps1

# macOS / Linux
# source .venv/bin/activate

# Install rockmoc in editable mode with plots + CLI + test tools
pip install -e ".[plot,cli,dev]"
```

### Minimal install (library only)

```bash
pip install -e .
```

### Optional chemistry stacks

```bash
pip install -e ".[cea]"        # RocketCEA (NASA CEA wrapper)
pip install -e ".[cantera]"
pip install -e ".[coolprop]"
pip install -e ".[all]"        # everything
```

### Check that it works

```bash
py -3.12 -c "import rockmoc; print(rockmoc.__version__)"
py -3.12 examples/quickstart.py
```

You should see MLN / ideal / TIC summaries and a CSV written in the current directory.

---

## 3. Mental model: how a design works

Every conventional design follows the same pattern:

```text
  gas model (γ, R, …)
        +
  DesignSpec (Rt, Me or ε or pe/pc, mesh density, …)
        │
        ▼
  designer function  →  DesignResult
        │
        ├── contour   (x, r, θ, M along the wall)
        ├── metrics   (ε, L, Cf_1d, angles, …)
        ├── net       (optional MOC field for plots)
        └── summary / export / plot
```

**Units (always SI in the API):**

| Quantity | Unit |
|----------|------|
| Length, radius | meters |
| Angles in code | **radians** (degrees appear in printed summaries) |
| Pressure | pascals (when used) |
| Temperature | kelvin (when used) |
| R | J/(kg·K) |

**Throat origin:** wall starts at \(x = 0\), \(r = R_t\).

---

## 4. Your first nozzle (step by step)

### Step A — Pick a gas

For a first run, constant γ is enough:

```python
from rockmoc import PerfectGas

# LOX/RP-1-ish engineering guess
gas = PerfectGas(gamma=1.25, R=350.0)

# Or air-like defaults
# gas = PerfectGas()  # γ=1.4, R≈287
```

**Tip:** `R` only matters for \(c^*\), \(I_{sp}\), and dimensional thermo. Contour *shape* for perfect-gas MOC is driven mainly by **γ** and the expansion target.

### Step B — Write a DesignSpec

Provide **exactly one** expansion target:

```python
from rockmoc import DesignSpec

spec = DesignSpec(
    gas=gas,
    throat_radius=0.05,       # 50 mm throat
    exit_mach=3.0,            # ← pick ONE of: exit_mach | area_ratio | pressure_ratio
    n_characteristics=48,     # 24–40 = quick; 80–160 = production quality
)
```

### Step C — Call a designer

```python
from rockmoc import design_mln

result = design_mln(spec, store_net=True)  # store_net=True if you want MOC plots
```

### Step D — Inspect, export, plot

```python
print(result.summary())                 # dict of key numbers
print(result.metrics.area_ratio)
print(result.metrics.length / spec.throat_radius)  # L / Rt

result.export_csv("my_nozzle.csv")
result.export_json("my_nozzle.json")

# Needs: pip install -e ".[plot]"
result.plot()                           # dark dashboard (contour + MOC + profiles)
# result.plot(show=False, save_dir="figures")
```

---

## 5. Choosing inputs: Me vs area ratio vs pe/pc

`DesignSpec` requires **exactly one** of:

| Field | Meaning | When to use |
|-------|---------|-------------|
| `exit_mach` | Design exit Mach \(M_e\) | You know the target Mach (most tutorials) |
| `area_ratio` | \(\varepsilon = A_e/A^*\) (axisym: \((R_e/R_t)^2\)) | You know geometric expansion from 1-D sizing |
| `pressure_ratio` | \(p_e / p_c\) | You size from chamber and exit static pressure |

```python
# All three of these are valid — pick one per DesignSpec

DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0)
DesignSpec(gas=gas, throat_radius=0.05, area_ratio=25.0)
DesignSpec(gas=gas, throat_radius=0.05, pressure_ratio=0.01)  # pe/pc
```

If you set two (or none), rockmoc raises `DesignError` with a clear message.

### Other useful DesignSpec fields

| Field | Default | Role |
|-------|---------|------|
| `n_characteristics` | 160 | Mesh density (fan rays). Higher → smoother, slower |
| `throat_mach` | 1.05 | Slightly supersonic start on the IVL |
| `throat_rc_over_rt` | 1.5 | Throat wall curvature / \(R_t\) (ideal / TIC arc) |
| `sharp_corner` | True | Sharp throat corner (classic MLN) |
| `truncate_fraction` | None | TIC only: fraction of ideal length (e.g. 0.8) |
| `delta` | 1.0 | 1 = axisymmetric, 0 = planar (advanced) |

---

## 6. Which nozzle type should I use?

```text
Need shortest high-ε vacuum nozzle?          →  design_mln
Need classical “full” contoured bell?        →  design_ideal
Need shorter than ideal, still contoured?    →  design_tic (truncate_fraction=0.7–0.9)
Need industry compact high-thrust shape?     →  design_rao (length_fraction≈0.8)
Need 2D / turbine-style MLN paper parity?    →  design_planar_mln_nasa
Need altitude-compensating plug?             →  design_spike / design_spike_truncated
Need “export this for review” package?       →  design_production
```

### Length intuition (same Me, γ, Rt)

```text
MLN  <  Rao TOP (~80% of 15° cone)  ≲  TIC (e.g. 80% ideal)  <  Ideal
```

Compare fairly:

```python
from rockmoc import compare_nozzle_designs, DesignSpec, PerfectGas

cmp = compare_nozzle_designs(
    DesignSpec(PerfectGas(1.25), throat_radius=0.05, exit_mach=3.0, n_characteristics=40)
)
print(cmp.summary_table())
```

---

## 7. Recipe book (copy-paste)

### 7.1 Axisymmetric MLN

```python
from rockmoc import PerfectGas, DesignSpec, design_mln

gas = PerfectGas(gamma=1.25, R=350.0)
spec = DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=64)
mln = design_mln(spec, store_net=True)
print(mln.summary())
mln.export_csv("mln.csv")
```

### 7.2 Ideal + TIC

```python
from rockmoc import PerfectGas, DesignSpec, design_ideal, design_tic

gas = PerfectGas(1.25, R=350.0)
base = dict(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=64)

ideal = design_ideal(DesignSpec(**base, throat_rc_over_rt=1.5))
tic = design_tic(DesignSpec(**base, truncate_fraction=0.80, throat_rc_over_rt=1.5))

print("ideal L/Rt", ideal.metrics.length / 0.05)
print("tic   L/Rt", tic.metrics.length / 0.05)
```

### 7.3 Rao TOP (thrust-optimized length)

```python
from rockmoc import PerfectGas, DesignSpec, design_rao

spec = DesignSpec(PerfectGas(1.25), throat_radius=0.05, exit_mach=3.5, n_characteristics=48)
rao = design_rao(spec, length_fraction=0.80)  # fraction of equivalent 15° cone length
print(rao.summary())
```

### 7.4 NASA planar MLN (TM-X-1502 sample)

```python
from rockmoc import design_planar_mln_nasa

# Paper sample case
r = design_planar_mln_nasa(
    exit_mach=4.05,
    gamma=1.36,
    dnu_deg=0.5,              # 0.1° = fine/slow; 0.5–1° = interactive
    throat_half_height=0.05,  # meters (planar half-height)
)
print(r.summary())
r.contour.export_csv("planar.csv")
```

### 7.5 Aerospike

```python
from rockmoc import design_spike, design_spike_truncated

spike = design_spike(
    Pc=5.0e6,          # chamber pressure [Pa]
    Tc=3200.0,         # chamber temperature [K]
    gamma=1.25,
    molar_mass=22.0,   # g/mol (Lee/Vernacchia convention)
    Re=0.15,           # shroud lip radius [m]
    expansion_ratio=25.0,
    Pa=0.0,            # ambient [Pa]; 0 = vacuum performance
    n_points=5000,
)
print(spike.summary())

# Practical truncated plug
trunc = design_spike_truncated(spike, length_fraction=0.25)
```

### 7.6 Area-ratio targeting (no Me)

```python
from rockmoc import PerfectGas, DesignSpec, design_mln

spec = DesignSpec(
    gas=PerfectGas(1.2),
    throat_radius=0.04,
    area_ratio=40.0,          # ε = Ae/A*
    n_characteristics=80,
)
r = design_mln(spec)
print("Design Me =", r.exit_mach_design)
print("Geometric ε =", r.metrics.area_ratio)
```

---

## 8. Reading the results

### `DesignResult` (most designers)

| Attribute | What it is |
|-----------|------------|
| `result.contour` | Wall polyline: `.x`, `.r`, `.theta`, `.M` |
| `result.metrics` | `area_ratio`, `length`, `exit_radius`, `Cf_vac_1d`, `max_wall_angle`, … |
| `result.exit_mach_design` | Design Me used |
| `result.net` | MOC field (if `store_net=True`) |
| `result.kind` | `"mln"`, `"ideal"`, `"tic"`, `"rao"`, … |
| `result.summary()` | Flat dict for printing / JSON |
| `result.export_csv(path)` | Wall table |
| `result.plot(...)` | Dark dashboard |

### Contour columns (CSV)

Typical columns: axial station `x` [m], radius `r` [m], wall angle, edge Mach (when available).

### Sanity checks after every design

```python
Rt = spec.throat_radius
assert result.contour.r[0] == Rt or abs(result.contour.r[0] - Rt) < 1e-9
assert result.contour.r[-1] > Rt                 # diverging
assert result.metrics.length > 0
assert abs(result.contour.theta[-1]) < 0.05       # ~axial exit for MLN/ideal-like
print("ε =", result.metrics.area_ratio)
print("L/Rt =", result.metrics.length / Rt)
print("θ_max deg =", result.summary()["max_wall_angle_deg"])
```

---

## 9. Student helpers (reports, gases, tables, cone)

### Fast homework path

```python
from rockmoc import DesignSpec, design_mln, design_conical, compare_nozzle_designs, example_gas

# Throat in millimeters, n=40 by default
spec = DesignSpec.student(Me=3.0, gas_name="lox_rp1", Rt_mm=50)
# or: DesignSpec.student(Me=3.0, gamma=1.25, Rt_mm=50)

mln = design_mln(spec)
print(mln.report(length_unit="mm", angle_unit="deg"))  # unit-aware summary
print(mln.explain())                                     # why this shape?
mln.lab_report("my_lab.md", student_name="Alex", course="AE-4XX")

cone = design_conical(spec, half_angle_deg=15.0)
cmp = compare_nozzle_designs(spec, include_conical=True)
print(cmp.to_markdown())          # paste into report
cmp.to_csv("compare.csv")
cmp.to_latex()                    # optional LaTeX tabular
```

### Example gases

```python
from rockmoc import example_gas, describe_example_gases, list_example_gases

print(describe_example_gases())
gas = example_gas("air")          # also: lox_rp1, lox_lh2, lox_ch4, …
```

These are **frozen teaching estimates**, not flight thermo.

### Known limits (always quote in labs)

```python
from rockmoc import KNOWN_LIMITS
for line in KNOWN_LIMITS:
    print("-", line)
```

Every `result.summary()` also includes a short `limitations` string.

### Notebooks & homework set

| Resource | Path |
|----------|------|
| Tutorials | `examples/notebooks/01_…` … `04_…` |
| Worked problems | [student_problems.md](student_problems.md) |
| Symbols | [glossary.md](glossary.md) |
| Mini explainers | [nozzle_explainers.md](nozzle_explainers.md) |

---

## 10. Plots and figures

### One-liner after design

```python
result = design_mln(spec, store_net=True)
result.plot()  # contour + MOC field + wall profiles + M/p
```

### Optional save (no GUI)

```python
result.plot(show=False, save_dir="figures", prefix="my_mln")
# → figures/my_mln_contour.png, …_moc_field.png, …
```

### Production dashboard

```python
prod = design_production(spec, kind="mln", apply_bl=True)
prod.plot(show=False, save_dir="figures", prefix="prod")
# also BL + altitude when available
```

### Standalone plot helpers

```python
from rockmoc.plot import (
    plot_nozzle_contour,
    plot_moc_field,
    plot_performance_vs_altitude,
    plot_of_ratio_sweep,
    of_ratio_sweep_rows,
    altitude_sweep_rows,
)

plot_nozzle_contour(result.contour, show=True, save_path=None)
plot_moc_field(result.contour, result.net, show=True)

rows = altitude_sweep_rows(
    result.contour, gas,
    throat_radius=0.05, exit_mach=result.exit_mach_design,
    pc=5e6, Tc=3200.0,
)
plot_performance_vs_altitude(rows, show=True)
```

### Generate the full demo gallery

```bash
py -3.12 tests/generate_all_graphs.py
# writes many PNGs under tests/graphs/  (gitignored by default)
```

---

## 11. Export for CAD / CFD

```python
# Wall only
result.export_csv("wall.csv")
result.export_json("wall.json")  # includes summary metadata

# Or from the contour object
result.contour.export_csv("wall_only.csv")
```

**Suggested workflow**

1. Design with rockmoc (start with `n_characteristics=40–80`).  
2. Mesh-converge (double `n`, compare \(L/R_t\) and wall shape).  
3. Export CSV → CAD (revolve wall) or CFD grid generator.  
4. Optional: `design_production(..., apply_bl=True)` and export **BL-corrected** wall.  
5. Final performance from CFD + real CEA thermo — rockmoc is the contour starter, not the last word.

---

## 12. Performance, altitude, and sweeps

### Integrated Cf / Isp (contour + 1-D edge)

```python
from rockmoc import integrated_performance

perf = integrated_performance(
    result.contour,
    gas,
    throat_radius=0.05,
    exit_mach=result.exit_mach_design,
    pc=5.0e6,
    Tc=3200.0,
    pa=0.0,            # vacuum
)
print(perf.summary())
# Cf_vac, Isp_vac, divergence_efficiency, …
```

### Altitude sweep

```python
from rockmoc import off_design_altitude_sweep

alt = off_design_altitude_sweep(
    result.contour, gas,
    throat_radius=0.05,
    exit_mach=result.exit_mach_design,
    pc=5e6, Tc=3200.0,
    altitudes_m=[0, 5e3, 10e3, 20e3, 40e3, 80e3],
)
for row in alt:
    print(row["altitude_m"], row["Cf"], row["Isp_s"])
```

### Field-integrated exit plane

```python
from rockmoc import field_integrated_cf

field = field_integrated_cf(
    result.contour, gas,
    throat_radius=0.05,
    exit_mach_design=result.exit_mach_design,
    pc=5e6, Tc=3200.0, pa=0.0,
    net=result.net,
)
print(field.summary())
```

### O/F or γ sweeps (for charts)

```python
from rockmoc.plot import of_ratio_sweep_rows, gamma_sweep_rows, plot_of_ratio_sweep

of_rows = of_ratio_sweep_rows(throat_radius=0.05, exit_mach=3.0, n_chars=24)
plot_of_ratio_sweep(of_rows, show=True, save_path="of_sweep.png")
```

---

## 13. Production pipeline (recommended for “final” contours)

One call stacks fine mesh, optional multi-γ, BL, field Cf, altitude, and certificates:

```python
from rockmoc import PerfectGas, DesignSpec, design_production

prod = design_production(
    DesignSpec(
        PerfectGas(1.25, R=350.0),
        throat_radius=0.05,
        exit_mach=3.0,
        n_characteristics=80,
    ),
    kind="mln",                         # mln | ideal | tic | rao
    apply_bl=True,
    multi_gamma=(1.25, 1.22, 1.18),    # optional γ_chamber, γ_throat, γ_exit
    pc=5.0e6,
    Tc=3200.0,
    n_chars_min=80,                     # raise for publication quality
    certificate=True,
)

print(prod.summary())

# Pieces you will use
wall_inviscid = prod.design.contour
wall_viscous  = prod.bl.corrected if prod.bl else wall_inviscid
print(prod.field.summary())
print("Altitude rows:", len(prod.altitude_sweep))

prod.plot(show=False, save_dir="prod_figures")
```

More detail: [production_pipeline.md](production_pipeline.md).

---

## 14. Thermochemistry (CEA / Cantera / CoolProp)

rockmoc MOC is still **perfect-gas / tabulated-γ** along the nozzle. Chemistry backends supply **better γ, R, Tc** at the design point (frozen surrogate — industry standard for preliminary MOC).

### RocketCEA example

```python
# pip install -e ".[cea]"
from rockmoc.thermo.rocketcea_backend import RocketCEAGas
from rockmoc import DesignSpec, design_mln

cea = RocketCEAGas(ox="LOX", fuel="RP1", of_ratio=2.3, Pc_psia=725.0)
gas = cea.solve()   # PerfectGas-like with chamber/throat properties

spec = DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.2, n_characteristics=64)
result = design_mln(spec)
```

If a backend is not installed, imports fail clearly — core designs never require CEA.

---

## 15. Command-line usage

```bash
pip install -e ".[cli]"

# MLN by exit Mach
rockmoc design --type mln --Me 3.0 --gamma 1.25 --Rt 0.05 -o mln.csv

# Ideal by area ratio
rockmoc design --type ideal --area-ratio 25 --Rt 0.04 -o ideal.csv

# TIC at 80% of ideal length
rockmoc design --type tic --Me 4.0 --truncate 0.8 --n 60 -o tic.csv

# Optional JSON summary
rockmoc design --type mln --Me 2.5 -o wall.csv --json wall_meta.json

rockmoc --version
```

| Flag | Meaning |
|------|---------|
| `-t / --type` | `mln` \| `ideal` \| `tic` |
| `-g / --gamma` | Ratio of specific heats |
| `--R` | Specific gas constant [J/(kg·K)] |
| `--Rt` | Throat radius [m] |
| `--Me` | Exit Mach |
| `-e / --area-ratio` | ε = Ae/A* |
| `--n` | Number of characteristics |
| `--truncate` | TIC fraction of ideal |
| `-o / --output` | Wall CSV path |

---

## 16. Mesh density: speed vs accuracy

| `n_characteristics` | Feel | Use when |
|--------------------:|------|----------|
| 16–24 | Seconds, coarse | Debugging scripts |
| 40–64 | Good interactive | Most exploration |
| 80–120 | Fine | Reports / CAD export |
| 160+ | Default production path | Final packages |

Planar NASA uses `dnu_deg` instead (degrees of Prandtl–Meyer step):

| `dnu_deg` | Notes |
|----------:|-------|
| 1.0 | Fast |
| 0.5 | Good balance |
| 0.1 | TM-X-1502 accuracy guidance (slow) |

**Rule of thumb:** double `n` once; if \(L/R_t\) and wall shape barely move, you are converged enough for preliminary design.

```python
from rockmoc import mesh_convergence_axisym, PerfectGas

rows = mesh_convergence_axisym(
    PerfectGas(1.25),
    throat_radius=0.05,
    exit_mach=3.0,
    n_list=[20, 40, 80],
)
for row in rows:
    print(row)
```

---

## 17. Validation: prove the install

### Fast unit tests

```bash
pip install -e ".[dev]"
pytest -q
```

### Deep literature campaign (writes error tables)

```bash
py -3.12 tests/run_deep_validation.py
```

Outputs (reference vs solver + **error %**):

```text
tests/data/validation_results/
  validation_results.json
  validation_results.csv
  validation_report.md
  sources.md
```

Open the CSV in Excel/Sheets to filter by regime (`thermo_isentropic`, `axisym_mln`, `rao_top`, …).

See [validation_deep.md](validation_deep.md) for sources and interpretation.

---

## 18. Common mistakes & FAQ

### “DesignError: Provide exactly one of exit_mach, area_ratio, or pressure_ratio”

You set zero or two targets. Leave the unused ones as `None`.

### “Design is slow”

Lower `n_characteristics` (try 32–48) or planar `dnu_deg=0.5`. Production mode forces a minimum mesh — pass `n_chars_min=48` for demos.

### “Plots crash with ImportError: matplotlib”

```bash
pip install -e ".[plot]"
```

### “I want degrees, not radians”

API stores radians. Convert when you need degrees:

```python
import math
theta_deg = math.degrees(result.contour.theta[-1])
```

Summaries already include degree fields where useful (`max_wall_angle_deg`, …).

### “ε from rockmoc ≠ my 1-D spreadsheet”

For axisym, ε is \((R_e/R_t)^2\). Check you did not compare to planar \(y_e/y_t\). Also confirm the same γ.

### “Can I use this for flight hardware?”

Use rockmoc to **generate and screen contours**, then:

1. Mesh-converge  
2. Attach real propellant thermo (CEA)  
3. CFD for shocks, separation, heat flux  
4. Structural / thermal design  

rockmoc is not a substitute for qualification testing.

### “Aerospike fails with expansion / ν errors”

Lee/Vernacchia needs \(\nu_e < \pi/2\) (finite PM turning). Extremely high ε with high γ can exceed that — lower ε or check γ.

### “CSV appears but git ignores it”

`.gitignore` blocks random `*.csv` exports so the repo stays clean. Your files still exist on disk; curated fixtures live under `tests/data/`.

### Windows: which Python?

```bash
py -3.12 -m pip install -e ".[plot,cli,dev]"
py -3.12 examples/quickstart.py
```

---

## 19. Where to go next

| Goal | Document / command |
|------|--------------------|
| Package overview | [README.md](../README.md) |
| Glossary & symbols | [glossary.md](glossary.md) |
| Homework problems | [student_problems.md](student_problems.md) |
| Nozzle explainers | [nozzle_explainers.md](nozzle_explainers.md) |
| Jupyter tutorials | `examples/notebooks/` |
| Production accuracy stack | [production_pipeline.md](production_pipeline.md) |
| First-principles / MOC equations | [first_principles.md](first_principles.md), [theory.md](theory.md) |
| NASA planar theory | [theory_nasa_tm_x_1502.md](theory_nasa_tm_x_1502.md) |
| Aerospike theory | [theory_aerospike_lee.md](theory_aerospike_lee.md) |
| Validation methodology | [validation.md](validation.md), [validation_deep.md](validation_deep.md) |
| Minimal script | `examples/quickstart.py` |
| Literature PDFs | `refs/README.md` |

### Minimal “happy path” checklist

- [ ] `pip install -e ".[plot,cli,dev]"`  
- [ ] `py -3.12 examples/quickstart.py`  
- [ ] Design one MLN with your Rt and Me  
- [ ] `result.plot()` and `result.export_csv(...)`  
- [ ] Optional: `design_production(...)` for BL + field Cf  
- [ ] `pytest -q` and/or `tests/run_deep_validation.py`  

---

*Happy designing. If something is confusing, open an issue with your `DesignSpec`, rockmoc version (`rockmoc.__version__`), and the full traceback.*
