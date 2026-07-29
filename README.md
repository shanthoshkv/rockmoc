# rockmoc

**Industry-oriented Method of Characteristics (MOC) rocket nozzle design for Python.**

Design conventional bells (MLN, ideal, TIC, Rao TOP/TOC), **planar NASA sharp-throat MLNs**, and **annular aerospikes** from literature-faithful algorithms — with pluggable thermochemistry, boundary-layer correction, field-integrated performance, dark-themed plots, and a multi-regime validation campaign against NASA/NACA/Rao sources.

| | |
|--|--|
| **Version** | see `src/rockmoc/_version.py` |
| **Python** | ≥ 3.10 |
| **Core deps** | NumPy, SciPy |
| **License** | MIT |
| **Author** | [shanthoshkv](https://github.com/shanthoshkv) |
| **Repo** | https://github.com/shanthoshkv/rockmoc |

---

## Table of contents

- [**How to use (full guide)**](docs/how_to_use.md) ← start here if you are new
- [What rockmoc is](#what-rockmoc-is)
- [Nozzle types](#nozzle-types)
- [Capabilities](#capabilities)
- [Install](#install)
- [Quickstart](#quickstart)
- [Dark plots & post-solve dashboard](#dark-plots--post-solve-dashboard)
- [Validation](#validation)
- [API map](#api-map)
- [Assumptions & limits](#assumptions--limits)
- [Project layout](#project-layout)
- [CLI](#cli)
- [References](#references)
- [Contributing / push checklist](#contributing--push-checklist)
- [License](#license)

> **New to rockmoc?** Start here:  
> - **[docs/how_to_use.md](docs/how_to_use.md)** — install, first design, student helpers, FAQ  
> - **[examples/notebooks/](examples/notebooks/)** — Jupyter tutorials (MLN, compare, altitude, NASA planar)  
> - **[docs/student_problems.md](docs/student_problems.md)** — worked homework  
> - **[docs/glossary.md](docs/glossary.md)** — symbols & units

---

## What rockmoc is

rockmoc is a **library + CLI** for *preliminary-to-production-engineering* nozzle contour design using the **method of characteristics**, not a full CFD solver.

It is built for:

- Propulsion students and researchers who need **reproducible MOC contours**
- Rocket teams who want **MLN / TIC / Rao / aerospike** options in one package
- Pipelines that export walls to **CAD / CFD** and attach **CEA / Cantera / CoolProp** γ later

**Coordinate system:** SI meters; geometric throat at \(x = 0\), wall radius \(r = R_t\); API angles in **radians** (summaries often print degrees).

---

## Nozzle types

| Type | Function(s) | Geometry | Typical use |
|------|-------------|----------|-------------|
| **Planar MLN** | `design_planar_mln_nasa` | 2D sharp-throat minimum-length, uniform exit | Turbine nozzles, educational MOC, NASA TM-X-1502 parity |
| **Axisym MLN** | `design_mln` | Axisymmetric sharp-throat minimum-length bell | Short high-ε vacuum stages, upper stages |
| **Ideal** | `design_ideal` | Circular throat arc + elongated canceling region | Full-length “ideal” contoured bell reference |
| **TIC** | `design_tic` | Truncated ideal contour (`truncate_fraction`) | Flight trade: length vs performance |
| **Rao TOP** | `design_rao` | Thrust-optimized parabolic / engineering TOP | Industry-standard compact high-thrust bells |
| **Rao TOC** | `design_rao_toc`, `design_rao_control` | Length-referenced / control-surface TOC variants | Advanced TOP–TOC studies |
| **Conical** | `design_conical` | Constant half-angle wall (default **15°**) | Textbook baseline vs MOC / TOP |
| **Aerospike (full)** | `design_spike` | Annular plug (Lee / Vernacchia) | Altitude-compensating engines |
| **Aerospike (truncated)** | `design_spike_truncated`, `truncated_spike_with_base` | Truncated spike + base pressure models | Practical plug nozzles |
| **Spike + cowl** | `design_spike_with_cowl`, `design_primary_cowl` | Primary cowl + spike stack | Integrated aerospike packages |
| **Production package** | `design_production` | Any of MLN/ideal/TIC/Rao + BL + field Cf + UQ | Recommended “ship the design” path |

### Contour family (same Me, γ)

Typical length ordering at fixed design point:

```text
MLN  <  Rao TOP (~0.8 × 15° cone)  ≲  TIC (e.g. 80% ideal)  <  Ideal
(+ conical 15° as length/angle baseline)
```

Use `compare_nozzle_designs(spec)` (includes cone by default) then:

```python
cmp = compare_nozzle_designs(spec)
print(cmp.to_markdown())   # lab report
cmp.to_csv("compare.csv")
```

**Student helpers:** `DesignSpec.student`, `example_gas`, `result.report()`, `result.lab_report()`, `explainer("mln")`, `KNOWN_LIMITS` — see [docs/how_to_use.md](docs/how_to_use.md#9-student-helpers-reports-gases-tables-cone).

---

## Capabilities

### Design & mesh

- **Axisymmetric MOC** unit processes (interior / axis / wall) with geometric source terms (Zucrow–Hoffman style)
- **Planar NASA MOC** (TM-X-1502): θ_max = ν_e/2, y_e = A/A*, mesh Δν control
- Dense defaults for accuracy (`n_characteristics` ~ 160, planar Δν down to 0.1°)
- Mesh studies: `mesh_convergence_planar`, `mesh_convergence_axisym`
- Optional **strict canceling wall** + robust schedule fallback
- **Isotropic area match** so geometric ε tracks 1-D A/A*(Me) without warping wall angles

### Thermochemistry

| Backend | Install | Role |
|---------|---------|------|
| `PerfectGas` | core | Constant-γ isentropic / PM |
| `VariableGammaGas` + expansion tables | core | Multi-γ along M (chamber → throat → exit) |
| RocketCEA | `pip install rockmoc[cea]` | Frozen chamber/throat γ, R, Tc from propellants |
| Cantera | `pip install rockmoc[cantera]` | Equilibrium / frozen mixtures |
| CoolProp | `pip install rockmoc[coolprop]` | Pure-fluid properties |

### Accuracy stack (production)

`design_production(...)` chains:

1. Fine MOC design (MLN / ideal / TIC / Rao)  
2. Expansion table (perfect-gas or multi-γ)  
3. Sauer-style **throat IVL** metadata  
4. Optional **BL δ\*** wall growth (TM-X-2343 spirit)  
5. **Field-integrated Cf / Isp** on the exit station  
6. Altitude off-design sweep  
7. Optional **mesh certificate / UQ** hooks  

### Performance & analysis

- 1-D ideal \(C_F\), \(I_{sp}\) (Sutton form)
- Exit-plane **divergence / kinetic** efficiencies (θ_e-based — no false wall-history penalty)
- `off_design_altitude_sweep` with ISA-style \(p_a(h)\)
- OF / γ / ε design **sweeps** for trade plots
- Monte Carlo / one-at-a-time sensitivity / mesh certificate

### Plotting (optional `rockmoc[plot]`)

Dark aerospace theme:

- Nozzle contour (mirrored walls, metrics box, wall-angle colormap)
- **MOC field**: C+, C−, expansion fan, canceling C−, wall nodes
- Wall profiles r(x), θ(x), M(x)
- Mach & p/p_c along wall
- Performance vs altitude
- O/F, γ, ε sweeps
- Family overlay + comparison bars
- Aerospike contour
- BL inviscid vs corrected

```python
result = design_mln(spec, store_net=True)
result.plot()                              # display dashboard
result.plot(show=False, save_dir="out")    # optional PNG save
```

### Export

```python
result.export_csv("wall.csv")
result.export_json("wall.json")
# Contour also: result.contour.export_csv / export_json
```

---

## Install

### From PyPI (when published)

```bash
pip install rockmoc
pip install "rockmoc[plot,cli]"   # figures + CLI
pip install "rockmoc[all]"        # all optional extras
```

### From GitHub

```bash
pip install "git+https://github.com/shanthoshkv/rockmoc.git"
```

### From a clone (development)

```bash
git clone https://github.com/shanthoshkv/rockmoc.git
cd rockmoc
pip install -e ".[plot,cli,dev]"

# Optional chemistry stacks
pip install -e ".[cea]"
pip install -e ".[cantera]"
pip install -e ".[coolprop]"
```

**Requirements:** Python ≥ 3.10, NumPy ≥ 1.24, SciPy ≥ 1.10.

Publishing notes for maintainers: see [PUBLISH.md](PUBLISH.md).

---

## Quickstart

For a slower, more detailed tour (inputs, mesh, exports, FAQ), see **[docs/how_to_use.md](docs/how_to_use.md)**.

### Axisymmetric MLN

```python
from rockmoc import PerfectGas, DesignSpec, design_mln

gas = PerfectGas(gamma=1.25, R=350.0)
spec = DesignSpec(
    gas=gas,
    throat_radius=0.05,      # m
    exit_mach=3.0,           # or area_ratio=... or pressure_ratio=...
    n_characteristics=80,
)
result = design_mln(spec, store_net=True)
print(result.summary())
result.export_csv("mln_wall.csv")
result.plot(show=True)       # needs rockmoc[plot]
```

### Ideal, TIC, Rao

```python
from rockmoc import design_ideal, design_tic, design_rao, DesignSpec, PerfectGas

spec = DesignSpec(PerfectGas(1.25), throat_radius=0.05, exit_mach=3.0, n_characteristics=64)

ideal = design_ideal(spec)
tic   = design_tic(DesignSpec(PerfectGas(1.25), throat_radius=0.05, exit_mach=3.0,
                              n_characteristics=64, truncate_fraction=0.80))
rao   = design_rao(spec, length_fraction=0.80)
```

### NASA planar MLN (TM-X-1502 sample)

```python
from rockmoc import design_planar_mln_nasa

# Paper sample: Me = 4.05, γ = 1.36
r = design_planar_mln_nasa(4.05, 1.36, dnu_deg=0.1, throat_half_height=0.05)
print(r.summary())
r.contour.export_csv("planar_mln.csv")
```

### Aerospike

```python
from rockmoc import design_spike, design_spike_truncated

spike = design_spike(
    Pc=5e6, Tc=3200.0, gamma=1.25, molar_mass=22.0,
    Re=0.15, expansion_ratio=25.0, Pa=0.0, n_points=5000,
)
print(spike.summary())
trunc = design_spike_truncated(spike, length_fraction=0.25)
```

### Production pipeline (recommended)

```python
from rockmoc import PerfectGas, DesignSpec, design_production

prod = design_production(
    DesignSpec(PerfectGas(1.25), throat_radius=0.05, exit_mach=3.0),
    kind="mln",                      # mln | ideal | tic | rao
    apply_bl=True,
    multi_gamma=(1.25, 1.22, 1.18), # optional γ_c, γ_t, γ_e
    pc=5e6,
    Tc=3200.0,
)
print(prod.summary())
prod.plot(show=False, save_dir="figures")
# prod.design.contour  — inviscid wall
# prod.bl.corrected    — BL-grown wall
# prod.field           — field Cf / Isp
```

### Fair comparison

```python
from rockmoc import DesignSpec, PerfectGas, compare_nozzle_designs

cmp = compare_nozzle_designs(
    DesignSpec(PerfectGas(1.25), throat_radius=0.05, exit_mach=3.0, n_characteristics=40)
)
print(cmp.summary_table())
```

### Chemistry-backed gas (optional)

```python
from rockmoc.thermo.rocketcea_backend import RocketCEAGas  # rockmoc[cea]

cea = RocketCEAGas(ox="LOX", fuel="LH2", of_ratio=5.0, Pc_psia=500.0)
gas = cea.solve()  # PerfectGas-like with frozen throat γ, R
```

---

## Dark plots & post-solve dashboard

```bash
pip install -e ".[plot]"
py -3.12 tests/generate_all_graphs.py   # writes tests/graphs/*.png
```

| Figure family | Contents |
|--------------|----------|
| Contour | Mirrored wall, throat/exit markers, metric box |
| MOC field | C+/C−, fan, canceling C−, wall & axis nodes |
| Wall profiles | r(x), θ(x), M(x) |
| Thermo edge | M(x) and p/p_c(x) |
| Altitude | C_F and I_sp vs h |
| Sweeps | O/F, γ, ε |
| Family | Overlay + bar comparison |
| Aerospike | Spike + lip radius |
| BL | Inviscid vs δ*-corrected wall |

Saving is always optional (`save_dir=None` → display only).

---

## Validation

rockmoc ships an offline, literature-anchored validation campaign.

### Run it

```bash
# Full unit + regression suite
pytest -q

# Deep multi-regime literature campaign → JSON / CSV / Markdown
py -3.12 tests/run_deep_validation.py
```

### Stored validation artifacts

| Path | Role |
|------|------|
| `tests/data/validation_catalog.json` | Reference targets, sources, tolerances |
| `tests/data/validation_results/validation_results.json` | Full ref vs solver |
| `tests/data/validation_results/validation_results.csv` | Flat table with **error_percent** |
| `tests/data/validation_results/validation_report.md` | Human-readable report |
| `tests/data/validation_results/sources.md` | Source index |
| `docs/validation.md` / `docs/validation_deep.md` | Methodology |

### Campaign snapshot (representative)

| Regime | Cases | What is checked | Typical error |
|--------|------:|-----------------|---------------|
| Isentropic / PM tables (γ=1.4) | 6 | ν, A/A*, pe/pc at M=1.5…5 | **≪ 0.02%** |
| Multi-γ round-trip | 5 | M→ν→M, M→A/A*→M | ~0 (machine) |
| NASA TM-X-1502 sample | 1 | θ_max=ν_e/2, y_e=A/A*, throat (0,1) | identity-level |
| NACA planar MLN suite | 5 | Me=2…4 design rules | ≪ 0.1% on rules |
| Axisym MLN multi-regime | 6 | ε vs 1-D, θ_e≈0, L/R_t band, diverging | ε match; θ_max ~4% vs planar half-ν by design |
| Rao TOP | 3 | L/L_15° cone, θ_n>θ_e, ε | few % on engineering L |
| Aerospike Lee | 2 | δ=π/2−ν_e, Me(ε), tip→0 | ≪ 1% identities |
| 1-D C_F (Sutton) | 5 | Closed-form C_F | ~0 |
| Family ordering | 1 | ideal > MLN; TIC < ideal; Rao compact | structural |

**Latest deep run:** **43/43 cases**, **170 metrics**, **12 cited sources**, median error ~0%, thermo table errors typically &lt; 0.02% (max ~4% only on documented axisym θ_max vs planar ν_e/2).

### Primary literature anchors (12)

1. **Anderson** — *Modern Compressible Flow* / γ=1.4 PM & isentropic tables  
2. **NACA Report 1135** (Ames, 1953) — compressible flow equations & tables  
3. **NASA TM-X-1502** — Vanco & Goldman (1968), planar sharp-edged MLN  
4. **NACA RM E8J12** — Shames & Seashore (1948), 2D sharp-edge design data  
5. **NASA TM-X-2343** — Goldman & Vanco (1971), MLN + BL displacement  
6. **Rao (1958)** — optimum thrust contour / TOP vs 15° cone  
7. **Lee / Vernacchia** — annular aerospike (δ = π/2 − ν_e)  
8. **Sutton & Biblarz** — ideal 1-D thrust coefficient  
9. **Zucrow & Hoffman** — axisymmetric MOC unit processes  
10. **Huzel & Huang** — conical nozzle practice (15° baseline)  
11. **Sauer / Kliegel** — transonic throat IVL  
12. **Shapiro** — Prandtl–Meyer function form & ν_max  

Catalog + CSV: `tests/data/validation_catalog.json`, `tests/data/validation_results/`.  
Local PDFs: `refs/*.pdf`.

### Also covered by pytest

- Planar MOC interior Riemann invariants  
- Physics regressions (wall angle vs geometry, canceling schedule, η_n = 1 at θ_e = 0)  
- Aerospike, production pipeline, thermo backends (optional skips)  
- Graph generation is separate (`tests/generate_all_graphs.py`)

---

## API map

```text
rockmoc
├── DesignSpec / DesignResult / ProductionResult
├── PerfectGas / VariableGammaGas
├── design_mln | design_ideal | design_tic
├── design_rao | design_rao_toc | design_rao_control
├── design_production
├── design_planar_mln_nasa
├── design_spike | design_spike_truncated | design_spike_with_cowl
├── correct_wall_bl | sauer_throat_ivl
├── integrated_performance | field_integrated_cf | off_design_altitude_sweep
├── compare_nozzle_designs
├── monte_carlo_mln | mesh_certificate | sensitivity_one_at_a_time
├── run_all_validation | run_deep_validation | write_validation_outputs
└── plot.*  (optional matplotlib dark theme)
```

Import surface is re-exported from the top-level `rockmoc` package for the main designers and utilities.

---

## Assumptions & limits

| Path | Assumptions |
|------|-------------|
| Planar NASA | 2D, inviscid, isentropic, constant γ |
| Axisym bells | Axisymmetric MOC; inviscid isentropic; γ from PerfectGas / frozen table / CEA |
| Canceling wall | Classic MLN (θ,ν) schedule + streamline integration + isotropic ε match (engineering hybrid, not full inverse CFD) |
| Aerospike | Lee/Vernacchia constant-γ isentropic annular model; needs ν_e &lt; π/2 |
| BL module | Integral compressible estimate (not full Cohen–Reshotko ODE suite) |
| Performance | Exit-station / 1-D edge models — validate with CFD before hardware |
| Chemistry | Frozen-γ surrogates for MOC; not multi-phase reacting CFD |

**Always:** mesh-converge → export wall → CFD / thermal / structural before flight hardware.

---

## Project layout

```text
rockmoc/
├── src/rockmoc/
│   ├── thermo/          # Perfect gas, variable-γ, CEA/Cantera/CoolProp
│   ├── moc/             # Equations, unit processes, planar NASA, axisym HF
│   ├── designers/       # MLN, ideal, TIC, Rao, production
│   ├── aerospike/       # Lee spike, base/cowl
│   ├── analysis/        # BL, throat IVL, atmosphere, UQ, compare
│   ├── performance/     # Metrics, integrated & field Cf
│   ├── geometry/        # Contour + export
│   ├── plot/            # Dark-themed figures
│   ├── validate/        # Literature suite + deep benchmark
│   └── cli.py
├── tests/
│   ├── data/
│   │   ├── validation_catalog.json
│   │   └── validation_results/   # ref vs solver + error %
│   ├── run_deep_validation.py
│   ├── generate_all_graphs.py
│   └── test_*.py
├── docs/
│   ├── how_to_use.md    # Friendly end-to-end usage guide
│   └── …                # Theory, validation, production notes
├── refs/                # Primary PDFs (NASA/NACA/Rao/…)
├── examples/
├── pyproject.toml
└── README.md
```

---

## CLI

```bash
pip install -e ".[cli]"

rockmoc design --type mln   --Me 2.5 --gamma 1.25 --Rt 0.05 -o wall.csv
rockmoc design --type ideal --area-ratio 25 --Rt 0.04 -o ideal.csv
rockmoc design --type tic   --Me 4.0 --truncate 0.8 -o tic.csv
rockmoc design --type rao   --Me 3.0 --gamma 1.25 --Rt 0.05 -o rao.csv
```

---

## References

- Zucrow, M. J. & Hoffman, J. D. — *Gas Dynamics* (axisym MOC)  
- Anderson, J. D. — *Modern Compressible Flow* / isentropic–PM tables  
- Vanco, M. R. & Goldman, L. J. — NASA TM-X-1502 (1968)  
- Goldman, L. J. & Vanco, M. R. — NASA TM-X-2343 (1971)  
- Shames, H. & Seashore, F. L. — NACA RM E8J12 (1948)  
- Rao, G. V. R. — Exhaust nozzle contour for optimum thrust (1958)  
- Lee, C. C. / Vernacchia — annular aerospike contour methods  
- Sutton, G. P. & Biblarz, O. — *Rocket Propulsion Elements*  
- Kliegel / Sauer — transonic throat approximations  

See `refs/README.md` and `docs/` for mapping from papers → modules.

---

## Contributing / release checklist

```bash
# 1. Install
pip install -e ".[plot,cli,dev]"

# 2. Tests
pytest -q

# 3. Deep validation (writes tests/data/validation_results/)
python tests/run_deep_validation.py

# 4. Build installable artifacts
python -m build
twine check dist/*

# 5. Optional: regenerate demo graphs (gitignored)
python tests/generate_all_graphs.py
```

Publish to PyPI: [PUBLISH.md](PUBLISH.md).

`.gitignore` excludes venvs, caches, `website/`, generated PNGs, secrets — while keeping `src/`, `tests/data/`, `docs/`, `examples/`.

---

## License

MIT — see [LICENSE](LICENSE).
