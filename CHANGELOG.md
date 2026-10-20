# Changelog

## 0.7.0

Physics fixes change numerical results for axisymmetric designs. Planar
(NASA TM-X-1502) results are unchanged.

### Fixed — physics
- **C− compatibility sign**: along C− (slope θ−μ), d(θ+ν)/dx = **+**sinθ sinμ / (r cos(θ−μ)).
  Unit processes carried θ−ν on C− and θ+ν on C+ (swapped); the axis point was
  reached along C+ (upstream). Verified against the exact conical source flow.
- **Axisymmetric MLN** rebuilt: predictor–corrector with full source terms (no blending /
  clamping); θ_max solved so the last corner C− reaches the axis at Me (θ_max < ν_e/2,
  e.g. 11.7° vs the old 25.9° for Me = 3, γ = 1.4); corner region resolved by C+ rows from
  the initial line; transition region solved as a Goursat problem to the uniform exit C+;
  wall = streamline carrying the throat mass flow. Exit area follows from mass conservation.
  Fan start consistent with the throat Mach (ν₀ = ν(M₀)).
- **Ideal contour** is now a true MOC solution (characteristic lattice from a Sauer
  transonic initial line along the circular-arc throat), replacing MLN + grafted arc +
  hand-tuned 1.55× stretch.
- **Rao TOC** is now built by MOC (control surface DE with Rao's conditions, mass-matched,
  wall by Goursat solve). The previous version ignored its kernel and drew a cubic Hermite.
- **Field-integrated Cf** now integrates momentum and pressure on the exit plane of the MOC
  net (the previous exit "ray" was synthetic); normalised by the field mass flow.
- **Sauer IVL** (`sauer_throat_ivl`) uses Sauer's solution as published (Mendelson,
  NASA TM X-53084, Eq. 2.6); the made-up "Hall factor" was removed.
- TIC planar area ratio; second-order (quadratic) mass-flow wall placement.

### Added
- `rockmoc.export`: CSV, CAD point list, DXF (R12) and OpenFOAM axisymmetric-wedge
  `blockMeshDict` (verified with `blockMesh` / `checkMesh` on OpenFOAM v2412).
- TOML design specs (`rockmoc run spec.toml`); CLI rewritten on `argparse` (Typer dropped).
- `analysis.grid_convergence` (Richardson extrapolation, observed order, GCI).
- `DesignResult.with_boundary_layer()` (δ* wall offset).
- Plots: Mach field fill, C+/C− families from real connectivity (`CharNet.c_minus`), wall
  coloured by Mach, uniform-exit region; normalised contour overlay with wall-angle inset;
  wall p/T/M properties; grid-convergence plot; Plotly interactive field (`[viz]`);
  light/publication theme; plots no longer mutate global rcParams; SVG default output.
- `MOCFallbackWarning`: numerical fallbacks are reported in `meta["warnings"]`.

### Changed
- Float kernels (`rockmoc._kernels`, Numba-compiled with `[fast]`): test suite 337 s → ~40 s.
- Variable-γ gas uses monotone cubic (PCHIP) tables with fast scalar lookups.
- Lazy top-level imports: `import rockmoc` 1.6 s → 0.03 s.
- `moc/axisym_hf.py` split into `_common`, `kernel_mln`, `transition`, `ideal`, `rao`.
- Single gas protocol (`MOCGas` incl. `gamma_at_mach`, `prandtl_meyer_max`); `FlowPoint` frozen.
- Version single-sourced from `_version.py`; extras `[fast]`, `[viz]`; Hypothesis property tests;
  CI on Linux + Windows, Python 3.10/3.12/3.13, nightly slow suite.

### Removed
- `moc/kernel.py`, `moc/strict_wall.py`, `mesh.wall_point_from_previous`,
  `mesh.wall_point_canceling` (under-determined / unused).

### Known limits
- Hall's higher-order transonic series is not implemented (Sauer only; accurate for Rc/Rt ≳ 2).
  For small Rc / low Me the ideal solver falls back to a uniform near-sonic start (reported).
- Rao TOC: very short bells can fall outside the feasible range on the lattice; the nearest
  feasible length is used with a warning.
- Axisymmetric MLN wall angle rises a few degrees past the corner before turning back
  (converged, shock-free); planar ideal shows a ~1° first-order bump after the arc attachment.

## 0.6.1
- Website / documentation links; plot documentation; validation catalog tweaks.
