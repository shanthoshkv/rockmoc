# rockmoc vs OpenFOAM: validation findings (v2, 2026-11-21)

45 CFD cases, `rhoCentralFoamAxi` (OpenFOAM v2412), 2.5° wedge (planar: 2-D), 8 MPI ranks each. v2 adds:
- a refined throat and mesh;
- wall-resolved viscous cases (y+ 0.2–0.8 on the production bell, about 2 at the throat);
- time-averaged fields;
- rockmoc fixes from v1 (`CHANGELOG.md`, Unreleased).

Results: **291/333 gated checks pass**. v1 was 225/335. Every one of the 42 remaining fails is explained below. Full tables and figures are in `results/REPORT.md`; regenerate with `py -3.12 cfdval.py compare`. v1 outputs are kept in `results/raw_v1`, `results/REPORT_v1.md` and `cases_v1/`.

## Vacuum Isp, rockmoc vs CFD

Sign in this file: rockmoc / CFD − 1 (positive = rockmoc high). The HTML/PDF report in `validation/` uses CFD / rockmoc − 1.

| nozzle | v1 (rockmoc v0.9) | v2 (fixed) |
|---|---|---|
| Rao 60 / 80 / 100 % | +3.3 / −0.6 / −2.7 % | **+0.02 / +0.04 / +0.03 %** |
| Rao TOC 80 % | −2.2 % | **+0.10 %** |
| Rao control 80 % | −1.4 % | **+0.12 %** |
| Cone 15° / 10° | −3.8 / −1.65 % | **−0.24 / −0.17 %** (mesh-independent) |
| Ideal M3, TIC 60 / 80 % | ±0.25 % | 0.00 / +0.15 / +0.01 % |
| MLN (Me 2–5, γ 1.15–1.4, n 10–80), production MLN | −0.06 % | ≤ 0.04 % |
| Planar MLN M2.4 (TM-X-1502) | −0.1 %, exit std 0.063 | **−0.03 %, exit std 0.008** |
| Full nozzles with chamber (ideal, TIC, MLN, Rao Rc 0.625/1/2) | — | −0.07 to +0.03 % |
| Production bell with BL correction (viscous, y+ ≈ 0.5) | −0.78 % | **−0.13 %** |
| MLN, uncorrected viscous wall | no friction in ref (+1.1 %) | **−0.14 %** |

Mesh independence (coarse / medium / fine, √2 per direction):
- Isp changes ≤ 0.02 %.
- MLN Isp GCI is 0.002 %.

## rockmoc fixes made from CFD evidence

1. **Divergence loss counted twice** (`performance/integrated.py`): λ and ∫cos²θ were both applied. Now the method is chosen automatically:
   - uniform exit → 1;
   - straight exit → λ (source flow);
   - otherwise the MOC wall-pressure integral.
2. **Rao TOP angles and parabola** (`designers/rao.py`):
   - Fixed a sign error in the θe correlation.
   - The bell is now an exact quadratic Bézier, with a Hermite fallback.
   - Fixes short bells: the MOC wall integral captures the lip over-expansion.
3. **Rao TOC** (`rao_toc.py`): it was built by stretching a TOP contour in x, which broke θ and Rc. It is now generated at the target length.
4. **Ideal/TIC throat arc**: polyline kinks of 5.9° are now resampled as an exact arc.
5. **Separation** (`analysis/separation.py`):
   - The Schmucker criterion is applied to the MOC wall pressure, not the 1-D pressure.
   - The vacuum part of Cf is kept separate so losses don't hit the ambient term.
6. **Boundary-layer loss.** rockmoc's wall-shear friction is right (CFD 1.49 % of thrust, rockmoc 1.52 %), but charging all of it as Isp loss overstates the loss by about 0.7 %.
   - BL-corrected wall: loss = exit control volume, P(ρu²θ − p_e δ*)/(p_c A_t), plus p_a on the displacement annulus at ambient.
   - Uncorrected wall: the core sees the effective area ratio, plus the momentum deficit (`delta_cf_uncorrected_vac`).
7. **Planar MLN wall** (`moc/planar_nasa.py`):
   - It was integrated with heuristic steps and rescaled to ε, which made it 7.8 % too long.
   - Now the wall comes from C+ / streamline intersections. Raw ε is within 0.02 % of 1-D, and L/h_t is 8.089 (independent MOC: 8.087).
8. **Transonic start line** (`moc/transonic.py`): Hall's 2nd-order series was used for Rc/Rt ≥ 1.5. It is 6 % high on the axis at 1.5 and has a non-physical minimum near 2.
   - Kliegel–Levine is now used up to Rc/Rt 4.
   - Throat-plane Mach now matches CFD within 0.3–2.4 % (axis) and 0.1–0.7 % (wall) at Rc/Rt 0.625–2. v1 was 3–9 % off.
9. **CLI default `-n` 40 → 80.** MLN axis error converges first-order in n: Mach rms 3.8 / 2.0 / 1.0 / 0.64 % at n = 10 / 20 / 40 / 80. The API default (160) was already fine.

## Verified correct

| what | evidence |
|---|---|
| MLN contour and field (Me 2–5, γ 1.15–1.4) | exit M within 0.2 %, exit std ≤ 0.014; MOC field Mach rms 0.3–0.5 %, flow angle < 0.4° |
| Axial M, p, T, ρ, u (MLN n = 80, ideal, TIC) | all axis gates pass: Mach rms 0.64 %, max 1.7 % |
| Discharge coefficient | 8 nozzles with a throat, Rc/Rt 0.625–2: ≤ 0.06 % inviscid, −0.32 % viscous D-5377 |
| Cuffel 1969 | Cd: CFD 0.981, rockmoc 0.982, measured 0.985; throat axis M 0.804 vs 0.80 measured |
| TN D-5377 wall pressure | 20 taps, rms log error 0.047 (cooled), within 0.1 for all runs |
| Real gas (LOX/RP-1, frozen table) | Isp −0.04 % |
| Viscous wall resolution | y+ 0.22–0.80 (median 0.46) on the production bell; CFD wall friction 1.49 % vs rockmoc 1.52 % |
| BL momentum thickness (production bell) | CFD/rockmoc 0.90 |
| Sea-level separation | Rao 80 % at 1 MPa: Cf +0.17 %; at 1.5 MPa: attached, as predicted (Cf −0.7 %). Cone 15° at 1 MPa: ε_sep 4.00 vs 3.90, Cf 1.135 vs 1.130 (−0.4 %, Kalt–Badal recovery correct). Schmucker p_sep/p_a 0.412 vs CFD 0.407 (Rao), 0.419 vs 0.448 (cone) |

## Remaining fails (42) and known offsets

- **MLN axis at n ≤ 40 (38 fails; includes frozen-gas M3 and the MLN mesh-convergence variants).**
  - This is rockmoc's own discretization (first order in n). It passes at n = 80, and the CLI default is now 80.
  - On the fine mesh the peak axis error grows (0.033 → 0.043) because CFD sharpens the kernel-end kink.
  - The gates are kept at n = 40 to document the trend.
- **Rao control 80 % (2 fails):** field Mach rms 1.4 %, flow angle rms 0.9° (v1: 1.8 %, 1.1°). Isp +0.12 %. Most likely a weak internal shock that the MOC net does not carry.
- **D-5377 cooled BL (2 fails):** CFD θ is 30–45 % below both the survey and rockmoc. rockmoc vs the survey is +5 to +25 %. This is a limit of the SST cooled-wall model in CFD, not a rockmoc defect.
- **Rao 80 % separation at 1 MPa: fixed by the Kliegel–Levine start line (fix 8).** The first v2 reference (Hall start) gave ε_sep 5.45 and Cf 1.041 (−4.2 %).
  - Now rockmoc predicts attached flow, with Cf 0.9955 vs CFD 0.9972 (+0.17 %).
  - The case is borderline: exit wall p/p_a is 0.377 against a Schmucker threshold of 0.3745 (0.6 % margin).
  - CFD separates only over the last 0.5 in ε (5.87 → 6.34), so thrust is unaffected.
  - The ε_sep check (−7.4 %) passes. With a flat bell lip, a 1 % wall-pressure change moves ε_sep by about 4 %.
- **Cone offset:** −0.2 % Isp, mesh-independent, because λ assumes source flow from a point. The MOC wall integral is worse (+0.6 % at 10°; at 15° the march covers only 20 % of the wall).

## CFD-side and design-assumption limits (informational checks)

- **MLN with a real converging section:**
  - Axis rms is 2.4 % against 1.0 % with a sonic start.
  - The MLN design assumes a uniform sonic line at a sharp corner; a real convergent gives a curved one.
  - Isp is unaffected (0.00 %).
- **Equilibrium-designed wall in frozen CFD:** −3.9 % Isp is a bound, not a match. The perfect-gas contour in frozen CFD gives +0.8 % Isp and misses Me by +1.9 %, so the variable-γ path matters.

## Not covered

- Aerospike (needs an external-flow mesher).
- Two-phase flow.
- Finite-rate or equilibrium chemistry inside CFD.
