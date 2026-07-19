# Validation status (v0.3)

## Automated tests

```text
45 passed, 1 skipped  (optional backend if not installed)
```

See generated **`docs/validation_report.md`** from `rockmoc.validate.validation_report_markdown()`.

| Suite | Status |
|-------|--------|
| Perfect-gas PM / area-Mach tables | Pass |
| Planar MOC invariants (δ=0) | Pass |
| NASA TM-X-1502 sample Me=4.05, γ=1.36 | Pass (θ_max=ν_e/2, y_e=A/A*) |
| NASA Δν mesh study | Pass |
| Aerospike Lee/Vernacchia identities | Pass |
| Axisym HF MLN / ideal / TIC | Pass |
| Axisym mesh study | Pass |
| Offline validation suite (`run_all_validation`) | Pass |
| RocketCEA / Cantera / CoolProp | Skip unless installed |

## Literature anchors

1. **NASA TM-X-1502** (NTRS 19680005278) — planar sharp-throat MLN  
2. **Lee 1963 / Vernacchia Spike Contour Algorithm** — annular aerospike  
3. **Anderson** compressible-flow tables — PM, A/A* (γ=1.4)  
4. **Zucrow & Hoffman** — axisymmetric MOC unit processes  
5. **Majdalani & Maickie** — Stodola inverse cross-check  

## Optional further work

- Digitize NACA RM E8J12 / Edelman contour *x* tables for pointwise planar regression  
- Bit-level float compare vs `mvernacc/nozzle_solver.py`  
- CEA Isp parity when RocketCEA is installed in CI  
- CFD smoke tests (out of core library)

## Honesty

rockmoc is a **traceable engineering design library**. Axisymmetric geometric ε is a MOC **output** (not forced equal to 1-D by default). Planar NASA sets y_e = A/A* per TM-X-1502 uniform-exit theory. Not a substitute for viscous CFD or hot-fire.
