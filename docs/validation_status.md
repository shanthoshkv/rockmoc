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

## Axisymmetric MOC corrections

- C− compatibility sign fixed: along C− (slope θ−μ), d(θ+ν)/dx = **+**sinθ sinμ / (r cos(θ−μ));
  invariants were also swapped in the unit processes. Verified against the exact conical source flow.
- Axisym MLN: θ_max solved so the last corner C− reaches the axis at Me (θ_max < ν_e/2, e.g.
  11.8° vs 24.9° planar for Me = 3, γ = 1.4); corner region and transition region (Goursat
  problem to the uniform exit C+) are solved, wall = streamline carrying the throat mass flow.
  Kernel mass conservation → 0.998 πRt² at N = 64; exit wall angle → 0 with N.
- Ideal contour: real MOC lattice from a Hall (2nd-order) transonic initial line along a circular-arc throat
  (uniform near-sonic start as fallback for small Rc / low Me), attachment chosen so the axis
  reaches Me, same transition + mass-flow wall. Exit ε reflects the discharge coefficient.
- Known limits: axisym MLN wall angle rises a few degrees past the corner before turning back
  (converged, shock-free); planar ideal shows a ~1° first-order bump after the arc attachment.

## Honesty

rockmoc is a **traceable engineering design library**. Axisymmetric geometric ε is a MOC **output** (not forced equal to 1-D by default). Planar NASA sets y_e = A/A* per TM-X-1502 uniform-exit theory. Not a substitute for viscous CFD or hot-fire.
