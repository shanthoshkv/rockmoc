# Changelog

## 1.0.1 — aerospike CFD tier, README

- New OpenFOAM tier T5 (`cfd_validation/spike.py`): Lee/Vernacchia annular spike at design, `pa = 0.3 pe`, `pa = 3 pe` and 25 % truncated. Cf within 0.7 % of CFD at design; overexpanded Cf under-predicted by 6.8 % (no altitude-compensation recompression in the fixed wall-pressure model).
- README rewritten for PyPI and GitHub (absolute links so images render on PyPI).

## 1.0.0 — fixes from the OpenFOAM validation campaign (`cfd_validation/`)

### Fixed
- **Divergence loss counted twice** (`performance.integrated`). `η_div · η_ke` modelled the
  same exit-angle loss twice. Against CFD, cones were 1.6–3.8 % low on Isp and bells 1.4–2.7 % low.
  The stack now applies one geometric factor `η_geom`. η_ke is still reported but no longer applied.
- **Rao TOP angle correlation inverted** (`designers.rao`). θ_e used to fall as the bell got
  shorter. It now rises, following the Huzel & Huang trends (60 % ≈ +5.5°, 100 % ≈ −3° vs 80 %).
  The parabola is now a true quadratic Bézier between the θ_n and θ_e tangents. The old cubic
  Hermite overshot θ_e by about 1°. Bells too short for a valid parabola fall back to the Hermite
  (`meta["parabola"]`).
- **Rao TOC stretched in x** (`designers.rao_toc`). The TOP contour was scaled axially after
  design, which left θ unscaled (stated θ_e 10.7°, real 5.6°) and made the Rc = 1.5 Rt throat arc
  an ellipse with Rc ≈ 5.6 Rt. The bell is now built at the target length (`design_rao_top(length_m=)`).
- **Ideal/TIC throat arc sampled by only 4–5 MOC nodes** (`designers.axisym_hf`). This put
  3–6° slope kinks into exported and meshed walls. The arc is now resampled with 40 exact points.
- **Separated sea-level thrust** (`analysis.separation`, `performance.integrated`):
  - Separation can be placed on the real wall pressure (`wall_eps`, `wall_p_pc`, `wall_M`).
  - `integrated_performance` applies losses to the vacuum part only, via
    `SeparationResult.Cf_vac_attached`.
  - Default recovery is Kalt–Badal (`separation_recovery=`).
  - CFD of Rao 80 % at p0/pa ≈ 10: v0.9 gave Cf 1.196, CFD gives 1.034.
- **BL-corrected nozzle loss** (`design_production(apply_bl=True)`): friction was charged as
  ∫τ_w dA (1.5 % Isp). The displaced wall still carries p_e, so the real loss is the exit
  control volume ΔCf = P(ρu²θ − p_e δ*)/(p_c A_t), plus p_a on the displacement annulus at
  ambient. New `BLCorrectionResult.meta["delta_cf_exit_vac"]`, `["exit_displacement_area_m2"]`;
  `exit_displacement_area=` on `field_integrated_cf` / `integrated_performance`.
  CFD at y+ ≈ 0.5: wall shear 1.49 % (rockmoc 1.52 %), net Isp 276.45 s; v0.9 274.31 s, now 276.09 s.
- **Uncorrected viscous wall**: new `BLCorrectionResult.meta["delta_cf_uncorrected_vac"]`, the
  Isp-equivalent loss when the wall is *not* BL-corrected (core sees (A_e − Pδ*_e)/(A_t − Pδ*_t),
  plus the momentum deficit). Pass it as `delta_cf_friction`. CFD 275.55 s; this model 275.17 s;
  ∫τ_w dA gave 274.31 s.
- **Planar MLN (`design_planar_mln_nasa`)**: the wall was integrated with heuristic arc steps and
  rescaled to ε, which made it 7.8 % too long at Me 2.4 (CFD exit Mach 2.45 on the axis, 2.20 at the wall).
  The wall now comes from C+ / streamline intersections of a standard unit-process net. The unforced
  MOC ε is within 0.02 % of 1-D, and L/h_t = 8.089 (independent MOC: 8.087). New `meta["area_ratio_moc_raw"]`.
- **Transonic start line**: Kliegel–Levine is now used up to Rc/Rt = 4 (it was 1.5). Throat-plane Mach
  from OpenFOAM (Rc/Rt = 0.625, 1, 1.5, 2) follows KL within 0.02 on the axis and wall. Hall's 2nd-order
  series is 6 % high on the axis at 1.5 (0.934 vs CFD 0.880) and has a non-physical minimum near Rc/Rt = 2.
  The ideal/TIC builders retry the series start at axis M 1.02 → 1.002 before falling back to a uniform start.
- CLI `-n` default 40 → 80 (TOML template 60 → 80). CFD: MLN axis Mach rms error converges first
  order in n (n = 10/20/40/80: 3.8/2.0/1.0/0.64 %); n ≥ 80 needed for < 1 %. API default stays 160.

### Added
- `performance.integrated.moc_thrust_efficiency`: inviscid thrust efficiency of any given wall,
  from a direct MOC march (`moc.analysis.analyze_wall`) plus a start-line and wall-pressure
  momentum balance. Results are cached.
- `integrated_performance(geometric="auto"|"divergence"|"moc", rc_over_rt=)`. `"auto"` selects:
  - θ_e < 0.5° → 1.
  - Straight exit → (1+cos θ_e)/2 (exact for source flow).
  - Otherwise the MOC wall integral, falling back to the divergence factor.
  - CFD (v1 meshes): cones, TIC, TOC and control-surface Rao are within 0.3 % Isp.

## 0.9.0 — real-gas flow, sweeps, units and export

### Added
- **Real-gas nozzle flow** (`thermo.cea_table.build_cea_table`, `--gas-model`, `[gas] model`):
  equilibrium or frozen-at-throat expansion table straight from RocketCEA, so
  γ_s = ρa²/p, T, ρ, a, molecular weight and the generalized Prandtl–Meyer angle
  ν = ∫√(M²−1) d ln V all vary along the nozzle (a constant-γ gas cannot do this). It is the
  default model whenever `[propellants]` is given; `perfect` keeps the old constant-γ gas and
  `multi_gamma` (TOML) blends γ_c/γ_t/γ_e. `--eps-max` bounds the table (default 400); designs
  beyond it are rejected, not clipped.
  - Validated: Mach, Isp (≤ 0.1 %) and c* (1e-5) reproduce RocketCEA's own values; an independent
    Cantera (h2o2) isentropic expansion of LOX/LH2 agrees on T_e (≤ 0.5 %) and V_e (≤ 0.15 %).
  - At the same ε = 25 (LOX/RP-1, 70 bar, O/F 2.4) a constant throat γ gives Me 3.656, L/Rt 23.3,
    Isp_vac 350.1 s; equilibrium 3.785 / 24.0 / 346.4 s (CEA: 346.4); frozen 4.042 / 25.5 / 332.8 s.
  - Boundary layer, field-integrated Cf and the Rao exit angle use the real-gas edge state
    (`thermo.edge_state`) when the gas has a table; constant-γ results are unchanged.
  - `VariableGammaGas`: p and T interpolated as cubic ln p / ln T (was linear: ~0.2 % Isp error);
    pe/pc → M inverts the gas's own p(M).
- **`rockmoc sweep`** / `rockmoc.sweep.run_sweep`: Cartesian grid over Me, area ratio, pe/pc, pe,
  pc, Tc, pa, altitude, O/F, gas model, γ, R, Rt, n, Rc/Rt, type, length, truncate (pressures
  accept `_bar/_mpa/_kpa/_psia`); parallel `--jobs`; failing cases land in an `error` column.
  Rows carry geometry, Cf, Isp, thrust, mdot, `gas_model`, `gamma_throat`, `gamma_exit`.
  `[sweep]` tables in a TOML spec run through `rockmoc run`. Sweep `gas_model=frozen,equilibrium`
  to bracket real-engine performance.
- **Performance from a design**: `DesignResult.performance(pc, Tc, pa)` → Cf, Isp, thrust, mdot;
  `[performance]` / `[propellants]` TOML tables; `--pc/--Tc/--pa/--altitude`;
  `thermo.propellant.of_props` (RocketCEA → γ, R, Tc).
- **Units**: `m, cm, mm, in, ft` for input (`DesignSpec.from_units`, `--in-units`, TOML
  `[design] units`) and output (`--units`); the solver stays SI.
- **Export**: `export_contour(path, n=, units=)` (`.csv`/`.txt`, chosen point count), `n` on
  `export_csv` / `export_cad_points`, `sldcrv` and `shell` (inner + outer offset) formats,
  `blockMeshDict` `units` (writes `scale`), `--stdout`, `--points`, `--thickness`;
  `Contour.offset`, `Contour.from_csv`; JSON export carries a units note.
- **CLI**: `rockmoc validate` (check a spec without solving), `rockmoc compare a.csv b.csv`
  (`analysis.compare.compare_contours`), `--pe-pc`, `--gas-model`, `--propellants/--of`.
- `separated_cf` and `psep_over_pa` exported from `rockmoc` / `rockmoc.analysis`; the separation
  docstrings now state the criteria are empirical screening correlations (±20–30 % scatter).
- Clear ImportError on Python < 3.10; PyPI publishing via trusted publishing (OIDC).
- Tests: `test_cea_table`, `test_sweep`, `test_moc_corner`, `test_separation`,
  `test_transonic_euler_conv` (+ units/export cases in `test_v07_features`).

### Changed
- With `[propellants]` / `--propellants` the default gas model is now **equilibrium real gas**
  (it was one constant γ and R). Pass `--gas-model perfect` for the old behaviour. Without
  RocketCEA it falls back to a constant-γ stand-in with a warning.

### Fixed
- **Euler throat solver now converges** (residual ~1e-7, was stalled at ~3e-2): the subsonic
  inflow BC copied interior static pressure (unstable at CFL 0.5); it is now characteristic
  (outgoing Riemann invariant). Cuffel / Kliegel–Levine results stay within published
  tolerances (Cuffel Rc/Rt = 0.625: axis M 0.825, Cd 0.981); ~2.5× faster.
- RocketCEA `get_Tcomb` returns °R; it was used as K in `plot/sweeps.py`,
  `thermo/expansion_table.py` and `thermo/rocketcea_backend.py` (Tc overstated ×1.8).

### Investigated, no change needed
- **Axisymmetric MLN "corner overshoot"** (wall angle ~15.5° vs corner θ_max ≈ 11.7°, γ = 1.4,
  Me = 3) is physics, not a numerical bug: the axisymmetric sources keep turning the wall outward
  after the corner. The source terms match the exact conical source flow (second order), and the
  peak is mesh-converged (15.45 / 15.48 / 15.50 / 15.50° at n = 16 / 32 / 64 / 128).
  `theta_max_rad` is the corner angle; `max_wall_angle_rad` is the peak. The "3° exit kink" does
  not reproduce (0.20 → 0.008° from n = 16 to 128).

### Known limits
- Equilibrium is the ideal upper bound and frozen the lower bound; finite-rate chemistry is not in
  the sweep path (see `thermo.kinetics`).
- The boundary-layer integral still uses the throat γ in its internal equations; only its edge
  states are real-gas.

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
