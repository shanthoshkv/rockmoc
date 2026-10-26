# Deep multi-regime validation

## How to run

```bash
py -3.12 tests/run_deep_validation.py
# or
py -3.12 -m rockmoc.validate.benchmark --out tests/data/validation_results
```

## Outputs (stored)

| File | Content |
|------|---------|
| `tests/data/validation_catalog.json` | Literature targets + **12 cited sources** + tolerances |
| `tests/data/validation_results/validation_results.json` | Full ref vs solver + error % |
| `tests/data/validation_results/validation_results.csv` | Flat table for spreadsheets |
| `tests/data/validation_results/validation_report.md` | Human-readable report |
| `tests/data/validation_results/sources.md` | Source index |

## Primary sources (12 cited)

| ID | Document |
|----|----------|
| `anderson_compressible` | Anderson, *Modern Compressible Flow* — γ=1.4 PM & isentropic tables |
| `naca_report_1135` | Ames Research Staff, NACA Report 1135 (1953) — compressible flow tables |
| `nasa_tm_x_1502` | Vanco & Goldman, NASA TM-X-1502 (1968) — planar sharp-throat MLN |
| `naca_rm_e8j12` | Shames & Seashore, NACA RM E8J12 (1948) — 2D sharp-edge design data |
| `nasa_tm_x_2343` | Goldman & Vanco, NASA TM-X-2343 (1971) — MLN + BL displacement |
| `rao_1958` | G.V.R. Rao (1958) — optimum thrust / TOP vs 15° cone |
| `lee_vernacchia_spike` | Lee 1963 / Vernacchia — aerospike δ=π/2−ν_e |
| `sutton_biblarz` | Sutton & Biblarz — ideal 1-D Cf |
| `zucrow_hoffman` | Zucrow & Hoffman — axisymmetric MOC unit processes |
| `huzel_huang` | Huzel & Huang — conical nozzle practice (15° baseline) |
| `sauer_kliegel` | Sauer/Hall — transonic throat IVL |
| `shapiro_gas_dynamics` | Shapiro — Prandtl–Meyer form & ν_max |

Local PDFs: `refs/`.

## Regimes covered

1. **Thermo isentropic** — M = 1.5…5, γ=1.4 (ν, A/A*, pe/pc)  
2. **ν_max (M→∞)** — Shapiro closed form  
3. **Multi-γ round-trip** — γ ∈ {1.2 … 1.4}  
4. **NASA TM-X-1502 sample** — Me=4.05, γ=1.36  
5. **NACA planar MLN suite** — Me = 2…4  
6. **Axisym MLN** — multi (γ, Me) contracts  
7. **Conical baseline** — Huzel L & η_div analytics  
8. **Rao TOP** — length fractions 0.7–0.8  
9. **Aerospike** — Lee identities  
10. **1-D Cf** — Sutton  
11. **BL δ\*** — TM-X-2343 spirit  
12. **Throat IVL** — Sauer/Hall  
13. **Family ordering** — MLN / TIC / ideal / Rao  

## Interpreting error %

```
error_percent = 100 * |solver − reference| / |reference|
```

- Thermo table metrics should be ≪ 0.1%.
- Design-rule identities (θ_max = ν_e/2, ε = A/A*) should be ≪ 1%.
- Engineering correlations (Rao L/L_15°) allow a few percent by design.
- Axisym θ_max vs pure planar ν_e/2 allows up to ~8% (documented mild axisym corner bump).
