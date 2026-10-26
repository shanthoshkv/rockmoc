# Changelog

## 0.8.0 — accuracy upgrades and known-limitation fixes

### Added
- **Kliegel–Levine transonic series** (`moc.transonic.kl_uv`, `kl_line`, `discharge_coefficient`):
  Hall's 3rd order re-expanded in 1/(1+R); used for 0.4 ≤ Rc/Rt < 1.5 (Hall above). Cuffel 1969
  Rc/Rt = 0.625: axis M 0.79 (0.8), Cd 0.982 (0.985). Axial scaling fixed by irrotationality.
- **2-D Euler throat solver** (`moc.transonic_euler`, Richardson-extrapolated first-order Rusanov FV)
  for Rc/Rt < 0.4 / steep inlets; matches KL at Rc/Rt = 2 (Cd 0.9970 vs 0.9962) and Cuffel at 0.625.
- **Bartz two-thickness boundary layer**: energy-integral equation (n = 0), H from two-layer 1/7
  profiles (`two_layer_thicknesses`), chamber-inlet start (`theta0`, `phi0`, subsonic contours),
  K_ax relaminarization switch (TN D-4788). TN D-5377 cooled-inlet θ now inside the measured band
  (was −38 to −43 %); adiabatic δ* +14 % worst (was +27 %).
- **Separation** (`analysis.separation`): Summerfield, Kalt–Badal, Schilling, Schmucker, Stark,
  Morrisette–Goldberg; truncated-nozzle or Kalt–Badal-recovery Cf (matches RocketCEA to 1e-6);
  TP-1207 full-flow/end-effect flags; `integrated_performance(..., separation="schmucker")`.
- **Internal-shock detector** (`moc.analysis.analyze_wall`): direct MOC on a given wall, cell-Jacobian
  sign test; flags TOP parabolas and arc-throat cones, not ideal contours.
- **Finite-rate kinetics** (`thermo.kinetics.kinetic_nozzle`, Cantera): quasi-1-D reacting nozzle;
  zero-rate → frozen, fast-rate → shifting equilibrium; equilibrium Isp within 0.1 % of CEA.
- **Two-phase lag model** (`performance.two_phase.two_phase_nozzle`): gas–particle velocity/thermal lag.
- **Validation**: Zebbiche axisymmetric MLN table (L/r* within 0.03 %), plus rows for all of the above;
  74/74 gated checks.

### Fixed
- Circular-import-safe lazy load kept; `correct_wall_bl` friction loss now integrates the divergent
  wall only (convergent friction belongs to c*).

## 0.8.0 (earlier items) — accuracy upgrades

Accuracy upgrades; numerical results change for BL-corrected, multi-γ and Rao designs.

### Added
- **Momentum-integral boundary layer** (`correct_wall_bl(method="integral")`, default):
  compressible von Kármán equation with pressure-gradient and axisymmetric divergence
  terms; Eckert/Sommer–Short reference-temperature Cf (within 3 % of Van Driest II,
  M ≤ 5); Crocco–Busemann shape factor `H = H_i T_w/T_e + r(γ−1)/2 M²`. Falls back to
  the flat-plate correlation (`method="empirical"`) with a `MOCFallbackWarning` if the
  ODE fails; the two agree exactly at zero pressure gradient.
- `wall_temperature` (cooled walls), `viscosity`, `theta0` options; Bartz transport
  defaults `μ = 1.184e-7 𝓜^0.5 T^0.6`, `Pr = 4γ/(9γ−5)` replace air Sutherland.
- **Skin-friction thrust loss** `BLCorrectionResult.delta_cf_friction`, and a JANNAF-style
  loss stack in `integrated_performance` / `field_integrated_cf`
  (`delta_cf_friction`, `eta_chemical`, `eta_two_phase`); production pipeline applies it.
- `rocketcea_backend.kinetic_efficiency_bounds`: frozen-at-throat / shifting-equilibrium
  bounds on η_chem.
- **Hall 2nd-order transonic initial line** (axisymmetric default, Sauer for planar).

- **Validation against published data** (`rockmoc.validate.experimental`,
  `docs/validation_experimental.md`, slow test): NASA TN D-5377 BL/wall pressure, Cuffel 1969
  transonic, NASA TM-107285 Rao-code bell, SP-8120 Rao chart, Angelino plug; 30/30 gated pass, 11 known-limit rows.
- **Website rebuilt** (`website/index.html`, `rockmoc.html`): heat-tint theme, live contour switcher and
  validation charts from real library output (`website/make_site_data.py` injects the data).

### Fixed
- **Hall throat series** restricted to Rc/Rt ≥ 1.5 (diverges below; Cuffel 1969 data). Sauer used
  for tighter throats, including the Rao default 0.382.
- **Circular import**: `from rockmoc import monte_carlo_mln` in a fresh interpreter.
- **Rao TOC control surface**: θ was eliminated through `acos`, singular at θ = μ, which
  the optimum DE crosses near L/L15 ≈ 0.8, so feasible lengths raised "not reachable" on
  some lattices (e.g. n = 20, 30, 60). Now a 2×2 Newton on (θ, M); converges monotonically in N.
- **Multi-γ table** now integrates the thermally perfect isentrope (h, p, V, ρV, ν) instead
  of applying constant-γ closed forms at the local γ; `VariableGammaGas` c* and Cf use the
  table state.

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
