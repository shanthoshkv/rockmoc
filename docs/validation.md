# Validation

> **Deep multi-regime campaign** (literature catalog + ref vs solver + error %):  
> see [`validation_deep.md`](validation_deep.md) and run  
> `py -3.12 tests/run_deep_validation.py` → `tests/data/validation_results/`.

## Analytic checks (automated)

| Check | Expected |
|-------|----------|
| \(\nu(M=2,\gamma=1.4)\) | \(26.380^\circ\) ± 0.01° |
| \(A/A^*(M=2,\gamma=1.4)\) | ≈ 1.6875 |
| \(A/A^*(M=3,\gamma=1.4)\) | ≈ 4.235 |
| PM inverse round-trip | \(\lvert\Delta M\rvert < 10^{-8}\) |
| Planar interior invariants | \(\theta\pm\nu\) conserved to \(10^{-8}\) |

## Design smoke tests

- MLN / ideal / TIC produce finite, positive wall radii and length  
- TIC length < ideal length at `truncate_fraction < 1`  
- Area-ratio targeting recovers design Mach  
- CSV / JSON export writable  

## Geometric design contracts (enforced)

| Contract | Behavior |
|----------|----------|
| Exit area ratio | Wall exit radius set so \(\varepsilon = (R_e/R_t)^2 = A/A^*(M_e)\) (1D isentropic) |
| Throat | Contour starts at \(r = R_t\), \(x = 0\) |
| Length | \(L/R_t\) kept in a physical band; non-finite / runaway nets fall back to a scheduled MLN wall |
| TIC | Axial length = `truncate_fraction` × ideal length |

## Known limitations (v0.1)

1. **Throat model** — MLN uses a sharp corner with a small geometric offset; ideal uses a simple circular arc + area–Mach estimate on the arc, not a full transonic CFD throat.  
2. **Canceling region** — hybrid MOC + streamline/schedule construction; exit uniformity is approximate vs. textbook hand-crafted nets at extreme \(M_e\). Kernel axis Mach may under-shoot design \(M_e\); exit size is still forced by isentropic \(\varepsilon(M_e)\).  
3. **Constant γ** — real motors need CEA / frozen–equilibrium properties.  
4. **Inviscid** — no boundary-layer displacement; wall is inviscid design wall.  
5. **Axis singularity** — source terms are regularized near \(r=0\); refine `n_characteristics` for stiff cases.  
6. **Length fidelity** — absolute \(L/R_t\) is engineering-grade for preliminary design; validate against CFD or published contours before hardware.  

## Recommended engineering workflow

1. Design with `rockmoc` (MLN or TIC).  
2. Mesh convergence: double `n_characteristics`, compare wall L2.  
3. Export wall → CAD / CFD.  
4. Correct with BL / CFD and CEA thermochemistry for final performance.
