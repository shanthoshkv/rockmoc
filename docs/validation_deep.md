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
| `tests/data/validation_catalog.json` | Literature targets + tolerances by regime |
| `tests/data/validation_results/validation_results.json` | Full ref vs solver + error % |
| `tests/data/validation_results/validation_results.csv` | Flat table for spreadsheets |
| `tests/data/validation_results/validation_report.md` | Human-readable report |
| `tests/data/validation_results/sources.md` | Source index |

## Primary sources

| ID | Document |
|----|----------|
| `anderson_naca_isentropic` | Isentropic / PM tables (γ=1.4), Anderson / NACA R-1135 lineage |
| `nasa_tm_x_1502` | Vanco & Goldman, NASA TM-X-1502 (1968) — planar sharp-throat MLN |
| `naca_rm_e8j12` | Shames & Seashore, NACA RM E8J12 (1948) — 2D sharp-edge design data |
| `nasa_tm_x_2343` | Goldman & Vanco, NASA TM-X-2343 (1971) — MLN + BL displacement |
| `rao_1958` | G.V.R. Rao (1958) — optimum thrust contour / TOP length fractions |
| `lee_vernacchia_spike` | Lee 1963 / Vernacchia spike contour — aerospike δ=π/2−ν_e |
| `sutton_1d_cf` | Sutton & Biblarz — ideal 1-D Cf |
| `zucrow_hoffman` | Zucrow & Hoffman — axisymmetric MOC unit processes |

Local PDFs: `refs/`.

## Regimes covered

1. **Thermo isentropic** — M = 1.5…5, γ=1.4 (ν, A/A*, pe/pc)
2. **Multi-γ round-trip** — γ ∈ {1.2, 1.25, 1.3, 1.36, 1.4}
3. **NASA TM-X-1502 sample** — Me=4.05, γ=1.36 planar MLN rules
4. **NACA planar MLN suite** — Me = 2…4
5. **Axisym MLN** — multi (γ, Me) design contracts (ε, θ_max, θ_e, L/Rt)
6. **Rao TOP** — length fraction of 15° cone, θ_n > θ_e
7. **Aerospike** — Lee identities, tip collapse, Me(ε)
8. **1-D Cf** — Sutton formula, vacuum + ambient
9. **Family ordering** — MLN / TIC / ideal / Rao length ordering

## Interpreting error %

```
error_percent = 100 * |solver − reference| / |reference|
```

- Thermo table metrics should be ≪ 0.1%.
- Design-rule identities (θ_max = ν_e/2, ε = A/A*) should be ≪ 1%.
- Engineering correlations (Rao L/L_15°) allow a few percent by design.
