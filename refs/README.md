# Reference library (`refs/`)

Primary literature used to ground **rockmoc** first principles. Prefer these
over secondary blogs when changing physics.

| File | Document | Use in rockmoc |
|------|----------|----------------|
| `NASA-TM-X-1502.pdf` | Vanco & Goldman, planar sharp-throat MLN MOC (1968) | Planar designer, Δν mesh |
| `NASA-TM-X-2343.pdf` | Goldman & Vanco, MLN + **BL displacement** (1971) | `correct_wall_bl` geometry rule |
| `NACA-RM-E8J12.pdf` | Shames & Seashore, sharp-edge throat design data (1948) | Validation of planar contours |
| `Rao_1958_Exhaust_Nozzle_Optimum_Thrust.pdf` | G.V.R. Rao, optimum thrust contour (1958) | Rao TOP / TOC designers |
| `Murthy_1965_nozzle_design.pdf` | Murthy, nozzle design theory (NASA CR / 1965) | Axisym MOC / Rao equations |
| `Korte_1997_aerospike.pdf` | Korte, multidisciplinary aerospike (1997) | Aerospike design context |
| `Spike_Contour_Algorithm.pdf` | Vernacchia / Lee plug contour notes | Annular aerospike |
| `IJERT_MLN_MOC.pdf` | MLN MOC design case study | Cross-check MLN procedure |
| `NASA-TN-D-5377_nozzle_BL.pdf` | Boldman et al., turbulent BL in M 4.4 cooled conical nozzle (1969) | `validate/experimental.py`: θ, δ*, wall p, throat p |
| `Cuffel_Back_Massier_1969_transonic.pdf` | Cuffel, Back & Massier, Rc/Rt = 0.625 transonic data (AIAA J. 1969, © AIAA) | Hall divergence guard, Sauer check |
| `NASA-TM-107285_nozzle_performance.pdf` | Arrington et al., Rao-code bell wall points (1996) | Rao TOC contour validation |
| `Aspire_TOP_nozzle_Rao_chart.pdf` | Newlands, TOP note reproducing SP-8120 θ_N/θ_e chart (2017) | Rao angle validation |
| `Kliegel_Levine_1969.pdf` | Transonic throat (Kliegel & Levine) | Sauer/Hall throat IVL |
| `Seitzman_Nozzle_Geometries.pdf` | GaTech nozzle geometry notes | Rao / MOC overview |

## Download script notes

Some URLs (Lee 1963 NTRS, Auburn thesis) may 404; re-fetch from NTRS search if needed.

## First-principles modules

| Physics | Code |
|---------|------|
| C± slopes, axisym sources | `src/rockmoc/moc/equations.py` |
| Unit processes (variable-γ) | `src/rockmoc/moc/mesh.py` |
| Throat IVL | `src/rockmoc/analysis/throat.py` |
| BL δ* add-on | `src/rockmoc/analysis/boundary_layer.py` |
| Rao control-surface TOC | `src/rockmoc/designers/rao_control.py` |
| NACA-style gold data | `tests/data/naca_rm_e8j12_gold.json` |
| Theory narrative | `docs/first_principles.md` |
