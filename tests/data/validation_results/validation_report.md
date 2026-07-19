# rockmoc deep multi-regime validation

**Generated (UTC):** 2026-07-19T07:53:00.747160+00:00

## Summary

| Cases | Passed | Failed | Metrics | Max err % | Mean err % | Median err % |
|------:|-------:|-------:|--------:|----------:|-----------:|-------------:|
| 34 | 34 | 0 | 151 | 4 | 0.203689 | 0 |

## Sources

- **anderson_naca_isentropic:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **nasa_tm_x_1502:** Vanco & Goldman, NASA TM-X-1502 (1968), Computer program for design of two-dimensional supersonic nozzle with sharp-edged throat
- **naca_rm_e8j12:** Shames & Seashore, NACA RM E8J12 (1948), Design data for graphical construction of two-dimensional sharp-edge-throat supersonic nozzles
- **nasa_tm_x_2343:** Goldman & Vanco, NASA TM-X-2343 (1971), Computer program for design of two-dimensional supersonic nozzle with boundary-layer correction
- **rao_1958:** G.V.R. Rao, Exhaust Nozzle Contour for Optimum Thrust, Jet Propulsion / ARS (1958)
- **lee_vernacchia_spike:** C.C. Lee (1963) annular aerospike; Vernacchia Spike Contour Algorithm notes
- **sutton_1d_cf:** Sutton & Biblarz, Rocket Propulsion Elements — ideal 1-D thrust coefficient
- **zucrow_hoffman:** Zucrow & Hoffman, Gas Dynamics — axisymmetric MOC unit processes

## Results by case

### `thermo_g1.4_M1.5` — **PASS**

- **Regime:** thermo_isentropic
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** Isentropic/PM table γ=1.4 M=1.5

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| nu_deg | 11.905 | 11.905209 | 0.001754 | 0.05 | Y | deg |  |
| area_ratio | 1.1762 | 1.176167 | 0.002801 | 0.05 | Y | - |  |
| pe_pc | 0.2724 | 0.272403 | 0.001126 | 0.1 | Y | - |  |
| Cf_vac_1d | — | 1.329625 | — | — | Y | - | solver 1-D Cf (no external table digit) |

### `thermo_g1.4_M2.0` — **PASS**

- **Regime:** thermo_isentropic
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** Isentropic/PM table γ=1.4 M=2.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| nu_deg | 26.38 | 26.379761 | 0.000906697 | 0.05 | Y | deg |  |
| area_ratio | 1.6875 | 1.6875 | 1.31582e-14 | 0.05 | Y | - |  |
| pe_pc | 0.1278 | 0.127805 | 0.003541 | 0.1 | Y | - |  |
| Cf_vac_1d | — | 1.423423 | — | — | Y | - | solver 1-D Cf (no external table digit) |

### `thermo_g1.4_M2.5` — **PASS**

- **Regime:** thermo_isentropic
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** Isentropic/PM table γ=1.4 M=2.5

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| nu_deg | 39.124 | 39.123564 | 0.001115 | 0.05 | Y | deg |  |
| area_ratio | 2.6367 | 2.636719 | 0.000711116 | 0.05 | Y | - |  |
| pe_pc | 0.05853 | 0.058528 | 0.003992 | 0.1 | Y | - |  |
| Cf_vac_1d | — | 1.50463 | — | — | Y | - | solver 1-D Cf (no external table digit) |

### `thermo_g1.4_M3.0` — **PASS**

- **Regime:** thermo_isentropic
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** Isentropic/PM table γ=1.4 M=3.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| nu_deg | 49.757 | 49.757347 | 0.000696876 | 0.05 | Y | deg |  |
| area_ratio | 4.2346 | 4.234568 | 0.000758012 | 0.05 | Y | - |  |
| pe_pc | 0.02722 | 0.027224 | 0.013533 | 0.1 | Y | - |  |
| Cf_vac_1d | — | 1.567815 | — | — | Y | - | solver 1-D Cf (no external table digit) |

### `thermo_g1.4_M4.0` — **PASS**

- **Regime:** thermo_isentropic
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** Isentropic/PM table γ=1.4 M=4.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| nu_deg | 65.785 | 65.78482 | 0.000273926 | 0.05 | Y | deg |  |
| area_ratio | 10.7188 | 10.71875 | 0.00046647 | 0.05 | Y | - |  |
| pe_pc | 0.006586 | 0.006586 | 0.001326 | 0.1 | Y | - |  |
| Cf_vac_1d | — | 1.651914 | — | — | Y | - | solver 1-D Cf (no external table digit) |

### `thermo_g1.4_M5.0` — **PASS**

- **Regime:** thermo_isentropic
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** Isentropic/PM table γ=1.4 M=5.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| nu_deg | 76.92 | 76.920216 | 0.000280172 | 0.05 | Y | deg |  |
| area_ratio | 25 | 25 | 2.84217e-14 | 0.05 | Y | - |  |
| pe_pc | 0.00189 | 0.00189 | 0.002031 | 0.1 | Y | - |  |
| Cf_vac_1d | — | 1.701035 | — | — | Y | - | solver 1-D Cf (no external table digit) |

### `thermo_multi_g1.2_M3.0` — **PASS**

- **Regime:** thermo_multi_gamma
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** PM & A/A* round-trip γ=1.2 Me=3.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| pm_roundtrip_M | 3 | 3 | 3.78956e-12 | 0.0001 | Y | - | M → ν → M |
| area_roundtrip_M | 3 | 3 | 0 | 0.001 | Y | - | M → A/A* → M |

### `thermo_multi_g1.25_M3.0` — **PASS**

- **Regime:** thermo_multi_gamma
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** PM & A/A* round-trip γ=1.25 Me=3.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| pm_roundtrip_M | 3 | 3 | 3.78956e-12 | 0.0001 | Y | - | M → ν → M |
| area_roundtrip_M | 3 | 3 | 0 | 0.001 | Y | - | M → A/A* → M |

### `thermo_multi_g1.3_M2.5` — **PASS**

- **Regime:** thermo_multi_gamma
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** PM & A/A* round-trip γ=1.3 Me=2.5

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| pm_roundtrip_M | 2.5 | 2.5 | 4.54747e-12 | 0.0001 | Y | - | M → ν → M |
| area_roundtrip_M | 2.5 | 2.5 | 4.54747e-12 | 0.001 | Y | - | M → A/A* → M |

### `thermo_multi_g1.36_M4.05` — **PASS**

- **Regime:** thermo_multi_gamma
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** PM & A/A* round-trip γ=1.36 Me=4.05

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| pm_roundtrip_M | 4.05 | 4.05 | 3.92553e-12 | 0.0001 | Y | - | M → ν → M |
| area_roundtrip_M | 4.05 | 4.05 | 1.68864e-12 | 0.001 | Y | - | M → A/A* → M |

### `thermo_multi_g1.4_M2.0` — **PASS**

- **Regime:** thermo_multi_gamma
- **Source:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage
- **Description:** PM & A/A* round-trip γ=1.4 Me=2.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| pm_roundtrip_M | 2 | 2 | 0 | 0.0001 | Y | - | M → ν → M |
| area_roundtrip_M | 2 | 2 | 0 | 0.001 | Y | - | M → A/A* → M |

### `nasa_tm_x_1502_sample` — **PASS**

- **Regime:** planar_mln_nasa
- **Source:** Vanco & Goldman, NASA TM-X-1502 (1968), Computer program for design of two-dimensional supersonic nozzle with sharp-edged throat
- **Description:** NASA TM-X-1502 sample Me=4.05 γ=1.36

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| theta_max_deg | 34.964558 | 34.964558 | 0 | 0.01 | Y | deg | θ_max = ν_e/2 (TM-X-1502) |
| area_ratio_1d | 12.981567 | 12.981567 | 0 | 0.0001 | Y | - | y_e/ht = A/A* for planar |
| ye_over_ht | 12.981567 | 12.981567 | 0 | 0.0001 | Y | - |  |
| throat_x | 0 | 0 | 0 | — | Y | m |  |
| throat_y | 1 | 1 | 0 | — | Y | m |  |
| xe_over_ht | — | 38.430906 | — | — | Y | - | physical length (mesh-dependent); TM-X-1502 discusses Δν sensitivity |

### `naca_planar_M2.0_dnu1.0` — **PASS**

- **Regime:** planar_mln_naca
- **Source:** Shames & Seashore, NACA RM E8J12 (1948), Design data for graphical construction of two-dimensional sharp-edge-throat supersonic nozzles
- **Description:** NACA E8J12-rule planar MLN Me=2.0 γ=1.4

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| theta_max_deg | 13.18988 | 13.18988 | 0 | 0.05 | Y | deg |  |
| ye_over_ht | 1.6875 | 1.6875 | 0 | 0.01 | Y | - | planar exit half-height = A/A* |
| area_ratio | 1.6875 | 1.6875 | 0 | 0.0001 | Y | - |  |
| L_over_ht | — | 6.208729 | — | — | Y | - | solver length |

### `naca_planar_M2.4_dnu1.0` — **PASS**

- **Regime:** planar_mln_naca
- **Source:** Shames & Seashore, NACA RM E8J12 (1948), Design data for graphical construction of two-dimensional sharp-edge-throat supersonic nozzles
- **Description:** NACA E8J12-rule planar MLN Me=2.4 γ=1.4

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| theta_max_deg | 18.373266 | 18.373266 | 0 | 0.05 | Y | deg |  |
| ye_over_ht | 2.4031 | 2.4031 | 0 | 0.01 | Y | - | planar exit half-height = A/A* |
| area_ratio | 2.4031 | 2.4031 | 0 | 0.0001 | Y | - |  |
| L_over_ht | — | 8.953531 | — | — | Y | - | solver length |

### `naca_planar_M3.0_dnu0.5` — **PASS**

- **Regime:** planar_mln_naca
- **Source:** Shames & Seashore, NACA RM E8J12 (1948), Design data for graphical construction of two-dimensional sharp-edge-throat supersonic nozzles
- **Description:** NACA E8J12-rule planar MLN Me=3.0 γ=1.4

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| theta_max_deg | 24.878673 | 24.878673 | 0 | 0.05 | Y | deg |  |
| ye_over_ht | 4.234568 | 4.234568 | 0 | 0.01 | Y | - | planar exit half-height = A/A* |
| area_ratio | 4.234568 | 4.234568 | 0 | 0.0001 | Y | - |  |
| L_over_ht | — | 14.867142 | — | — | Y | - | solver length |

### `naca_planar_M3.5_dnu0.5` — **PASS**

- **Regime:** planar_mln_naca
- **Source:** Shames & Seashore, NACA RM E8J12 (1948), Design data for graphical construction of two-dimensional sharp-edge-throat supersonic nozzles
- **Description:** NACA E8J12-rule planar MLN Me=3.5 γ=1.4

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| theta_max_deg | 29.264879 | 29.264879 | 0 | 0.05 | Y | deg |  |
| ye_over_ht | 6.789621 | 6.789621 | 0 | 0.01 | Y | - | planar exit half-height = A/A* |
| area_ratio | 6.789621 | 6.789621 | 0 | 0.0001 | Y | - |  |
| L_over_ht | — | 22.439111 | — | — | Y | - | solver length |

### `naca_planar_M4.0_dnu0.5` — **PASS**

- **Regime:** planar_mln_naca
- **Source:** Shames & Seashore, NACA RM E8J12 (1948), Design data for graphical construction of two-dimensional sharp-edge-throat supersonic nozzles
- **Description:** NACA E8J12-rule planar MLN Me=4.0 γ=1.4

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| theta_max_deg | 32.89241 | 32.89241 | 0 | 0.05 | Y | deg |  |
| ye_over_ht | 10.71875 | 10.71875 | 0 | 0.01 | Y | - | planar exit half-height = A/A* |
| area_ratio | 10.71875 | 10.71875 | 0 | 0.0001 | Y | - |  |
| L_over_ht | — | 33.277264 | — | — | Y | - | solver length |

### `axisym_mln_g1.2_M2.5` — **PASS**

- **Regime:** axisym_mln
- **Source:** Zucrow & Hoffman, Gas Dynamics — axisymmetric MOC unit processes
- **Description:** Axisym MLN γ=1.2 Me=2.5 n=40

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| area_ratio | 3.420533 | 3.420533 | 1.47228e-11 | 0.05 | Y | - | ε_geom vs A/A*(Me) |
| theta_max_vs_half_nu_e_deg | 24.230755 | 25.199985 | 4 | 8 | Y | deg | corner ≈ ν_e/2 (axisym mild offset allowed) |
| theta_exit_deg | 0 | 0 | 0 | — | Y | deg | axial exit design |
| L_over_Rt | — | 3.251889 | 0 | — | Y | - | band [1.5, 25.0] |
| exit_radius_m | 0.092473 | 0.092473 | 7.3686e-12 | 0.05 | Y | m |  |
| Cf_vac_1d | — | 1.579038 | — | — | Y | - |  |
| diverging | 1 | 1 | 0 | — | Y |  | Re > Rt |

### `axisym_mln_g1.25_M3.0` — **PASS**

- **Regime:** axisym_mln
- **Source:** Zucrow & Hoffman, Gas Dynamics — axisymmetric MOC unit processes
- **Description:** Axisym MLN γ=1.25 Me=3.0 n=40

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| area_ratio | 5.831867 | 5.831867 | 1.42246e-11 | 0.05 | Y | - | ε_geom vs A/A*(Me) |
| theta_max_vs_half_nu_e_deg | 29.706395 | 30.894651 | 4 | 8 | Y | deg | corner ≈ ν_e/2 (axisym mild offset allowed) |
| theta_exit_deg | 0 | 0 | 0 | — | Y | deg | axial exit design |
| L_over_Rt | — | 4.344173 | 0 | — | Y | - | band [1.5, 25.0] |
| exit_radius_m | 0.120746 | 0.120746 | 7.11438e-12 | 0.05 | Y | m |  |
| Cf_vac_1d | — | 1.648728 | — | — | Y | - |  |
| diverging | 1 | 1 | 0 | — | Y |  | Re > Rt |

### `axisym_mln_g1.3_M3.5` — **PASS**

- **Regime:** axisym_mln
- **Source:** Zucrow & Hoffman, Gas Dynamics — axisymmetric MOC unit processes
- **Description:** Axisym MLN γ=1.3 Me=3.5 n=40

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| area_ratio | 9.109847 | 9.109847 | 0 | 0.05 | Y | - | ε_geom vs A/A*(Me) |
| theta_max_vs_half_nu_e_deg | 33.159042 | 34.485403 | 4 | 8 | Y | deg | corner ≈ ν_e/2 (axisym mild offset allowed) |
| theta_exit_deg | 0 | 0 | 0 | — | Y | deg | axial exit design |
| L_over_Rt | — | 5.485957 | 0 | — | Y | - | band [1.5, 25.0] |
| exit_radius_m | 0.150913 | 0.150913 | 0 | 0.05 | Y | m |  |
| Cf_vac_1d | — | 1.680032 | — | — | Y | - |  |
| diverging | 1 | 1 | 0 | — | Y |  | Re > Rt |

### `axisym_mln_g1.4_M2.0` — **PASS**

- **Regime:** axisym_mln
- **Source:** Zucrow & Hoffman, Gas Dynamics — axisymmetric MOC unit processes
- **Description:** Axisym MLN γ=1.4 Me=2.0 n=40

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| area_ratio | 1.6875 | 1.6875 | 1.31582e-14 | 0.05 | Y | - | ε_geom vs A/A*(Me) |
| theta_max_vs_half_nu_e_deg | 13.18988 | 13.717476 | 4 | 8 | Y | deg | corner ≈ ν_e/2 (axisym mild offset allowed) |
| theta_exit_deg | 0 | 0 | 0 | — | Y | deg | axial exit design |
| L_over_Rt | — | 2.161682 | 0 | — | Y | - | band [1.5, 25.0] |
| exit_radius_m | 0.064952 | 0.064952 | 0 | 0.05 | Y | m |  |
| Cf_vac_1d | — | 1.423423 | — | — | Y | - |  |
| diverging | 1 | 1 | 0 | — | Y |  | Re > Rt |

### `axisym_mln_g1.4_M3.0` — **PASS**

- **Regime:** axisym_mln
- **Source:** Zucrow & Hoffman, Gas Dynamics — axisymmetric MOC unit processes
- **Description:** Axisym MLN γ=1.4 Me=3.0 n=40

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| area_ratio | 4.234568 | 4.234568 | 1.07809e-11 | 0.05 | Y | - | ε_geom vs A/A*(Me) |
| theta_max_vs_half_nu_e_deg | 24.878673 | 25.87382 | 4 | 8 | Y | deg | corner ≈ ν_e/2 (axisym mild offset allowed) |
| theta_exit_deg | 0 | 0 | 0 | — | Y | deg | axial exit design |
| L_over_Rt | — | 3.93657 | 0 | — | Y | - | band [1.5, 25.0] |
| exit_radius_m | 0.10289 | 0.10289 | 5.39518e-12 | 0.05 | Y | m |  |
| Cf_vac_1d | — | 1.567815 | — | — | Y | - |  |
| diverging | 1 | 1 | 0 | — | Y |  | Re > Rt |

### `axisym_mln_g1.4_M4.0` — **PASS**

- **Regime:** axisym_mln
- **Source:** Zucrow & Hoffman, Gas Dynamics — axisymmetric MOC unit processes
- **Description:** Axisym MLN γ=1.4 Me=4.0 n=40

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| area_ratio | 10.71875 | 10.71875 | 1.01589e-11 | 0.05 | Y | - | ε_geom vs A/A*(Me) |
| theta_max_vs_half_nu_e_deg | 32.89241 | 34.208106 | 4 | 8 | Y | deg | corner ≈ ν_e/2 (axisym mild offset allowed) |
| theta_exit_deg | 0 | 0 | 0 | — | Y | deg | axial exit design |
| L_over_Rt | — | 6.23701 | 0 | — | Y | - | band [1.5, 25.0] |
| exit_radius_m | 0.163698 | 0.163698 | 5.08662e-12 | 0.05 | Y | m |  |
| Cf_vac_1d | — | 1.651914 | — | — | Y | - |  |
| diverging | 1 | 1 | 0 | — | Y |  | Re > Rt |

### `rao_top_g1.25_M3.0_lf0.8` — **PASS**

- **Regime:** rao_top
- **Source:** G.V.R. Rao, Exhaust Nozzle Contour for Optimum Thrust, Jet Propulsion / ARS (1958)
- **Description:** Rao TOP γ=1.25 Me=3.0 L/L_15°=0.8

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| L_over_L15cone | 0.8 | 0.8 | 0 | 5 | Y | - | Rao length fraction of 15° cone |
| area_ratio | 5.831867 | 5.831867 | 0 | 0.5 | Y | - |  |
| theta_n_deg | — | 25.81435 | — | — | Y | deg |  |
| theta_e_deg | — | 11.425525 | — | — | Y | deg |  |
| theta_n_gt_theta_e | 1 | 1 | 0 | — | Y |  | Rao attachment angle > exit angle |

### `rao_top_g1.4_M3.5_lf0.8` — **PASS**

- **Regime:** rao_top
- **Source:** G.V.R. Rao, Exhaust Nozzle Contour for Optimum Thrust, Jet Propulsion / ARS (1958)
- **Description:** Rao TOP γ=1.4 Me=3.5 L/L_15°=0.8

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| L_over_L15cone | 0.8 | 0.8 | 0 | 5 | Y | - | Rao length fraction of 15° cone |
| area_ratio | 6.789621 | 6.789621 | 0 | 0.5 | Y | - |  |
| theta_n_deg | — | 26.182857 | — | — | Y | deg |  |
| theta_e_deg | — | 11.154254 | — | — | Y | deg |  |
| theta_n_gt_theta_e | 1 | 1 | 0 | — | Y |  | Rao attachment angle > exit angle |

### `rao_top_g1.25_M4.0_lf0.7` — **PASS**

- **Regime:** rao_top
- **Source:** G.V.R. Rao, Exhaust Nozzle Contour for Optimum Thrust, Jet Propulsion / ARS (1958)
- **Description:** Rao TOP γ=1.25 Me=4.0 L/L_15°=0.7

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| L_over_L15cone | 0.7 | 0.7 | 0 | 5 | Y | - | Rao length fraction of 15° cone |
| area_ratio | 20.644259 | 20.644259 | 1.72092e-14 | 0.5 | Y | - |  |
| theta_n_deg | — | 30.110415 | — | — | Y | deg |  |
| theta_e_deg | — | 8.544183 | — | — | Y | deg |  |
| theta_n_gt_theta_e | 1 | 1 | 0 | — | Y |  | Rao attachment angle > exit angle |

### `aerospike_lee_0` — **PASS**

- **Regime:** aerospike
- **Source:** C.C. Lee (1963) annular aerospike; Vernacchia Spike Contour Algorithm notes
- **Description:** Lee spike ε=25.0 γ=1.25

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| delta_rad | 0.105253 | 0.105253 | 0 | 0.01 | Y | rad | δ = π/2 − ν_e |
| exit_mach | 4.153776 | 4.153776 | 1.2188e-12 | 1 | Y | - | Me from ε (Majdalani / area–Mach) |
| expansion_ratio | 25 | 25 | 0 | 0.01 | Y | - |  |
| R_tip_near_zero | 0 | 6.12128e-09 | inf | — | Y | m | full spike tip → axis |
| Isp_s | — | 306.705 | — | — | Y | s | solver performance |
| Cf | — | 1.799821 | — | — | Y | - |  |

### `aerospike_lee_1` — **PASS**

- **Regime:** aerospike
- **Source:** C.C. Lee (1963) annular aerospike; Vernacchia Spike Contour Algorithm notes
- **Description:** Lee spike ε=12.0 γ=1.3

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| delta_rad | 0.334138 | 0.334138 | 0 | 0.01 | Y | rad | δ = π/2 − ν_e |
| exit_mach | 3.744466 | 3.744466 | 6.54665e-12 | 1 | Y | - | Me from ε (Majdalani / area–Mach) |
| expansion_ratio | 12 | 12 | 0 | 0.01 | Y | - |  |
| R_tip_near_zero | 0 | 0 | 0 | — | Y | m | full spike tip → axis |
| Isp_s | — | 221.965 | — | — | Y | s | solver performance |
| Cf | — | 1.300588 | — | — | Y | - |  |

### `cf_1d_g1.4_M2.0_pa0.0` — **PASS**

- **Regime:** performance_1d
- **Source:** Sutton & Biblarz, Rocket Propulsion Elements — ideal 1-D thrust coefficient
- **Description:** 1-D Cf γ=1.4 Me=2.0 pa/pc=0.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| Cf | 1.423423 | 1.423423 | 0 | 0.05 | Y | - | Sutton ideal Cf |
| pe_pc | 0.127805 | 0.127805 | 0 | 1e-09 | Y | - |  |
| area_ratio | 1.6875 | 1.6875 | 0 | 1e-09 | Y | - |  |

### `cf_1d_g1.4_M3.0_pa0.0` — **PASS**

- **Regime:** performance_1d
- **Source:** Sutton & Biblarz, Rocket Propulsion Elements — ideal 1-D thrust coefficient
- **Description:** 1-D Cf γ=1.4 Me=3.0 pa/pc=0.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| Cf | 1.567815 | 1.567815 | 0 | 0.05 | Y | - | Sutton ideal Cf |
| pe_pc | 0.027224 | 0.027224 | 0 | 1e-09 | Y | - |  |
| area_ratio | 4.234568 | 4.234568 | 0 | 1e-09 | Y | - |  |

### `cf_1d_g1.25_M3.0_pa0.0` — **PASS**

- **Regime:** performance_1d
- **Source:** Sutton & Biblarz, Rocket Propulsion Elements — ideal 1-D thrust coefficient
- **Description:** 1-D Cf γ=1.25 Me=3.0 pa/pc=0.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| Cf | 1.648728 | 1.648728 | 0 | 0.05 | Y | - | Sutton ideal Cf |
| pe_pc | 0.023078 | 0.023078 | 0 | 1e-09 | Y | - |  |
| area_ratio | 5.831867 | 5.831867 | 0 | 1e-09 | Y | - |  |

### `cf_1d_g1.25_M3.0_pa0.01` — **PASS**

- **Regime:** performance_1d
- **Source:** Sutton & Biblarz, Rocket Propulsion Elements — ideal 1-D thrust coefficient
- **Description:** 1-D Cf γ=1.25 Me=3.0 pa/pc=0.01

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| Cf | 1.590409 | 1.590409 | 0 | 0.05 | Y | - | Sutton ideal Cf |
| pe_pc | 0.023078 | 0.023078 | 0 | 1e-09 | Y | - |  |
| area_ratio | 5.831867 | 5.831867 | 0 | 1e-09 | Y | - |  |

### `cf_1d_g1.2_M4.0_pa0.0` — **PASS**

- **Regime:** performance_1d
- **Source:** Sutton & Biblarz, Rocket Propulsion Elements — ideal 1-D thrust coefficient
- **Description:** 1-D Cf γ=1.2 Me=4.0 pa/pc=0.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| Cf | 1.85415 | 1.85415 | 0 | 0.05 | Y | - | Sutton ideal Cf |
| pe_pc | 0.003237 | 0.003237 | 0 | 1e-09 | Y | - |  |
| area_ratio | 28.355261 | 28.355261 | 0 | 1e-09 | Y | - |  |

### `family_ordering_mln_tic_ideal_rao` — **PASS**

- **Regime:** family_ordering
- **Source:** Zucrow & Hoffman, Gas Dynamics — axisymmetric MOC unit processes
- **Description:** Length ordering γ=1.25 Me=3.0

| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |
|--------|----------:|-------:|--------:|------:|:----:|------|------|
| L_mln_m | — | 0.211959 | — | — | Y | m |  |
| L_tic_m | — | 0.250361 | — | — | Y | m |  |
| L_ideal_m | — | 0.312951 | — | — | Y | m |  |
| L_rao_m | — | 0.211223 | — | — | Y | m |  |
| ideal_gt_mln | 1 | 1 | 0 | — | Y |  | ideal longer than MLN |
| tic_lt_ideal | 1 | 1 | 0 | — | Y |  | TIC shorter than ideal |
| rao_lt_ideal | 1 | 1 | 0 | — | Y |  | Rao TOP typically shorter than full ideal |
