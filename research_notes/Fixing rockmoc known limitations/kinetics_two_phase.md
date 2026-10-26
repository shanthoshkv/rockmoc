# Finite-rate chemistry (kinetic) and two-phase losses in rocket nozzles, for rockmoc

> Research caveat: most primary sources I tried (the TDK/ODK manuals on NTRS, DTIC AD-A004666 "Specific Impulse Losses in Solid Propellant Rockets", the Stanford AA284a lecture) are scanned or compressed PDFs that the fetch tool could not parse, or they returned 403. Primary-source text on exact ODK, Bray, Hoglund, SPP and Hermsen equations was **not** retrieved. Equations below are marked as one of:
> - **[SOURCED]**: seen in a fetched source.
> - **[DERIVED]**: derived here from the standard conservation laws. These are safe to implement, but check them against the TDK manual.
> - **[RECALLED, UNVERIFIED]**: from background knowledge. Do not hard-code these without checking the cited original.

## 1. One-dimensional kinetics (ODK): governing equations, Cantera implementation, throat singularity, integration strategy

### Takeaway
ODK is the JANNAF quasi-1D finite-rate nozzle module, a companion to ODE, TDK, TBL and TDE. Its physics is a steady, inviscid, quasi-1D reacting flow along A(x). You can implement it in Cantera by integrating u, T and Y_k with an implicit ODE solver. The frozen-sonic singularity at the throat can be avoided in either of two ways: assume equilibrium up to just past the throat and start the kinetic integration at a slightly supersonic point, or march in pressure through the transonic region. I found no ready-made Cantera kinetic-nozzle script.

### Cited Findings
- ODK is one of five programs that the JANNAF (ex-ICRPG) Liquid Rocket Performance Committee selected in 1967. The other four are ODE (1-D equilibrium), TDK (2-D kinetics), TBL (turbulent boundary layer) and TDE (2-D equilibrium) — [NTRS TDK search results / 19860019159 context](https://ntrs.nasa.gov/citations/19860019159)
- TDK uses "a Two-Dimensional Method of Characteristics solution with fully coupled finite rate kinetics for axially symmetric nozzles", with optional frozen or equilibrium solutions. LPP is a superset of TDK — [Dunn & Coats, NTRS 19960029260](https://ntrs.nasa.gov/citations/19960029260) (PDF: https://ntrs.nasa.gov/api/citations/19960029260/downloads/19960029260.pdf)
- The ODK and TDK modules share a linkage with the ODE module. In 1984 they were extended to electrons and ionic species, with reaction rates taken from Jensen & Jones, "Reaction Rate Coefficients for Flame Calculations", Combust. Flame 32:1–34 (1978) — [Nickerson, Improved TDK progress report, NASA CR-171107 (N86-28631), PDF](https://ntrs.nasa.gov/api/citations/19860019159/downloads/19860019159.pdf)
- That report says the SSME is "a very demanding test case" because of wall-induced shocks and its "near equilibrium thermodynamic state of the exhaust" — [same PDF](https://ntrs.nasa.gov/api/citations/19860019159/downloads/19860019159.pdf)
- Other TDK documents on NTRS: the Engineering and Programming Manual for the 12/73 TDK revision (Ultrasystems), a shock-capability final report (NAS8-34974, March 1984), and [19850011050](https://ntrs.nasa.gov/api/citations/19850011050/downloads/19850011050.pdf) and [19770078415](https://ntrs.nasa.gov/api/citations/19770078415/downloads/19770078415.pdf). The CR number of the 1973 manual could not be read from the OCR — [NTRS search results](https://ntrs.nasa.gov/citations/19850011050)
- RocketIsp cites the TDK 1985 engineering manual, Dunn & Coats "Nozzle Performance Predictions Using the TDK 97 Code" (AIAA), *Liquid Rocket Thrust Chambers* ch. 17, and CPIA 246 (the JANNAF Rocket Engine Performance Prediction and Evaluation Manual). It says "JANNAF only controls a version of TDK from 1993" — [RocketIsp nozzle_eff](https://rocketisp.readthedocs.io/en/latest/nozzle_eff.html)
- NASA LSENS (Radhakrishnan & Bittker, RP-1328/1329/1330) solves steady 1-D inviscid reacting flow with **pressure or area assigned** as a function of distance or time. It uses LSODE as the stiff integrator and also gives sensitivity coefficients. It is a documented reference implementation for the same problem — [NTRS 19950008642](https://ntrs.nasa.gov/citations/19950008642)
- NATA (a Fortran code) solves quasi-1D nozzle flow for frozen, equilibrium and finite-rate nonequilibrium chemistry. It can serve as a physics reference — [CORE record](https://core.ac.uk/works/24923146)
- Cantera's example gallery has a "Converging-Diverging Nozzle" example and an isentropic-flow example under thermodynamics. It also has 1-D flame and reactor machinery. No official finite-rate nozzle example was found — [Cantera example tags](https://cantera.org/3.2/_tags/C.html), [Cantera 1-D docs](https://cantera.org/dev/reference/onedim/index.html)

### Inferences
**Governing equations [DERIVED].** These are steady, inviscid and adiabatic, with ideal-gas mixture, species k = 1..K, molecular weights W_k, molar production rates ω̇_k (Cantera `net_production_rates`, kmol/m³/s), and frozen sound speed a_f² = γ_f R T / W̄ with γ_f = cp/cv frozen.
```
continuity:  ρ u A = ṁ = const
momentum:    ρ u du/dx = −dp/dx
energy:      h(T,Y) + u²/2 = h0        (mass-specific h, includes formation)
species:     dY_k/dx = ω̇_k W_k / (ρ u)
state:       p = ρ R_u T / W̄
```
Eliminating dp, dρ and dT gives the working form, with M = u/a_f the frozen Mach number:
```
du/dx = u/(1 − M²) · [ −(1/A) dA/dx + Σ_k ( W̄/W_k − h_k/(c_p T) ) dY_k/dx ]
dT/dx = −( u du/dx + Σ_k h_k dY_k/dx ) / c_p          (h_k = mass-specific species enthalpy)
dp/dx = −ρ u du/dx
```
Here h_k/(c_p T) uses mass-specific h_k and c_p; `gas.partial_molar_enthalpies / W_k` gives h_k. The bracket is the "effective area change" caused by heat release and mole change. In the subsonic chamber the chemistry term is near zero.

**Sonic singularity.** At M = 1 the bracket must vanish. With finite-rate chemistry the bracket vanishes slightly downstream of the geometric throat, where dA/dx > 0. ṁ is an eigenvalue, so you cannot simply integrate in A(x) through M = 1. Implementable options, simplest first:
1. **Equilibrium-through-throat start (recommended for rockmoc).** Rocket chambers are near equilibrium up to the throat (high p and T). Use CEA or Cantera shifting equilibrium up to a station at A/A* ≈ 1.01–1.1 on the supersonic side, or at M ≈ 1.05–1.2. Set the composition and state there, then integrate the ODEs above in x with the supersonic A(x) from the MOC contour. This approach also gives ṁ and c*. I believe this is the classical ODK approach, but the TDK manual was unreadable to me, so it is unverified. A per-reaction or Bray-type check (section 2) can move the switch point upstream if the flow is not near equilibrium at the throat.
2. **Pressure as the independent variable.** In the transonic region, prescribe p and integrate du/dp = −1/(ρu) and dY_k/dp = ω̇_k W_k /(ρ u) · dx/dp. Then A = ṁ/(ρu) comes out as a result, and the minimum of A defines the throat. This requires prescribing dx/dp, i.e., a residence-time model; LSENS offers an "assigned pressure" mode for this. Use it to cross M = 1, then switch to assigned area.
3. **Shooting on ṁ.** Iterate ṁ (or the inlet Mach) so that the numerator and denominator vanish together. Then step across the singularity with L'Hôpital or linear extrapolation from M ≈ 0.98 to 1.02. This is the most expensive and fragile option.

**Cantera specifics [DERIVED].** Use `ct.Solution(mech)`. In the RHS, call `gas.TPY = T, p, Y`, then read `gas.net_production_rates`, `gas.molecular_weights`, `gas.cp_mass`, `gas.partial_molar_enthalpies`, `gas.mean_molecular_weight` and `gas.density`. Integrate the state vector [u, T, Y_1..K] with `scipy.integrate.solve_ivp(method='LSODA' or 'BDF')`, since the system is stiff. Alternatively, march the Lagrangian residence time with a `ct.IdealGasReactor` and update the state from the 1-D relations, but the direct ODE is simpler. Converge ṁ against rockmoc's existing CEA throat values.

**Isp from the solution.** Isp_vac = (u_e + p_e A_e/ṁ)/g0. Then η_kin = Isp_ODK/Isp_ODE (shifting), and RocketIsp's fraction is fracKin = (Isp_ODK − Isp_frozen)/(Isp_eq − Isp_frozen). The fraction slots into rockmoc's existing frozen/equilibrium bounds.

**Mechanisms (judgment, unverified).** For H2/O2, use a modern H2/O2 set (for example Burke 2012 or the H2/O2 subset of GRI-3.0). For CH4, GRI-Mech 3.0, which was validated at lower pressure. For RP-1, a surrogate or a reduced exhaust-recombination set: H+H+M, H+OH+M, O+O+M, O+H+M, CO+OH⇌CO2+H, H2+OH⇌H2O+H, O+H2⇌OH+H, H+O2⇌OH+O, OH+OH⇌H2O+O. In expansion, recombination reactions dominate. For N2O4/MMH, add N/NO chemistry (e.g., GRI NOx).

### Gaps
- The exact TDK/ODK equation set, the numerical method, the throat start logic and the default reaction/rate table (the JANNAF "standard" reaction set) were not retrieved. The NTRS PDFs are scans; get CR-178628 / TDK 1985 manual (Nickerson, Coats, Bartz) and CPIA 246 manually.
- No open-source Cantera quasi-1D kinetic nozzle script was found. The Cantera "Converging-Diverging Nozzle" example could not be confirmed as finite-rate; it is likely equilibrium or isentropic.

## 2. Bray sudden-freezing criterion

### Takeaway
Bray (1959, J. Fluid Mech.) showed that for a single dissociation–recombination reaction, nozzle flow leaves equilibrium over a short distance and then freezes almost immediately. Composition then stays roughly constant downstream. The freeze point is where the rate at which equilibrium composition changes along the flow equals the one-way recombination rate. The method works well for steep, single-reaction cases and degrades for slow (low-angle) expansions and multi-reaction rocket exhaust.

### Cited Findings
- Bray, "Atomic recombination in a hypersonic wind-tunnel nozzle", J. Fluid Mech. 6 (1959). Deviations from equilibrium occur over a very wide range of rate constants, and "once such a deviation has begun, the gas will very rapidly 'freeze'". The dissociation fraction then stays almost constant — [Cambridge abstract](https://resolve.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/atomic-recombination-in-a-hypersonic-windtunnel-nozzle/A180D2661B68DD9CB329906B6F33D23D)
- Bray's 1960 Southampton PhD thesis on the same topic includes approximate solutions and an optimum nozzle shape. It is open access — [PDF](https://eprints.soton.ac.uk/458375/1/38084.pdf)
- Watson (1961, Caltech): "suddenness" depends on both rate constants and geometry. Bray noted that "for nozzles with very small expansion angles, the freezing stretches out over quite some length and is far from sudden" — [Caltech authors record](https://authors.library.caltech.edu/records/rz0ga-7sp02)
- Porter's 1963 Caltech thesis: modified near-equilibrium criteria give results "substantially equivalent" to Bray's criterion — [Porter 1963 thesis PDF](https://thesis.library.caltech.edu/7035/3/Porter_jw_1963ocr.pdf)

### Inferences
- **Criterion [RECALLED, UNVERIFIED in exact form].** Take a reaction A2 + M ⇌ 2A + M with dissociation fraction α. Integrate equilibrium (α = α_e(x)) until
  `| dα_e/dt | = | u dα_e/dx |  ≈  R_r  = k_r ρ² α_e² (·)/W-scaling`,
  i.e., until the magnitude of the equilibrium composition's Lagrangian rate of change equals the one-way (recombination) rate evaluated at the equilibrium state. Downstream of that point, freeze the composition. Equivalently, freezing occurs where the ratio of flow time to chemical time is about 1. Check the exact definition (one-way rate vs. forward rate, and which constants) in Bray's thesis PDF above before coding.
- **Multi-species generalization.** A practical rockmoc approach: run shifting equilibrium along A(x) from rockmoc's contour and compute u(x). For each key recombination reaction, compare |u d[X]_e/dx| with the reaction's one-way rate from Cantera (`forward_rates_of_progress` / `reverse_rates_of_progress`). The freeze station is the first x where the dominant reactions fail the test. Isp is then computed by expanding frozen from that station: CEA/RocketCEA supports freezing at a given area ratio (`frozen=1, frozenAtThroat=0`, or the CEA "nfz" input point). This gives a cheap η_kin estimate between the current bounds.
- Accuracy relative to full kinetics: Bray is a single-reaction approximation. For H2/O2 and hydrocarbon exhaust with multiple coupled recombinations and heat release, expect it to bracket ODK rather than match it. Validate against the full ODE integration from section 1.

### Gaps
- The exact Bray equation and its quantitative accuracy against ODK for rocket propellants were not retrieved. The Southampton thesis PDF should contain the derivation.

## 3. Typical kinetic-loss magnitudes (η_kin) by propellant, chamber pressure and size

### Takeaway
Kinetic loss runs from about 0.1% to 10% of Isp. Thruster size (throat radius) is the dominant parameter, followed by chamber temperature, Pc, composition and ε. The loss is largest for small, low-Pc, high-ε engines. For LOX/RP-1 at moderate to high Pc it is very small. For LOX/LH2, CEA can exceed TDK by more than 3 s at high ε and low MR.

### Cited Findings
- "Kinetics loss is typically in the range of 0.1 to 10% of Isp." RocketIsp fits fracKin with MLP regressors on TDK results. Top parameters for the combined group are Rthrt, TcCham, Pc, HCL, *CO and eps, and "thruster size (Rthrt) is the most important factor" for most groups — [RocketIsp nozzle efficiencies](https://rocketisp.readthedocs.io/en/latest/nozzle_eff.html)
- Definitions: η_kin = Isp_with_loss/Isp_without, and η_kin' (fracKin) = the fraction of the gap between frozen and equilibrium, with 0 = frozen and 1 = equilibrium — [RocketIsp](https://rocketisp.readthedocs.io/en/latest/nozzle_eff.html)
- TDK is the JANNAF standard partly because "unlike CEA and RPA, the TDK code can calculate the kinetics loss", which matters "especially in small thrust, low pressure, high area ratio engines." LOX/LH2 shows close agreement in most cases, except "at high area ratio, low MR, where CEA gets more than 3 seconds higher than TDK." For LOX/RP-1 the chemical-kinetics impact is "very small" — [RocketCEA: CEA vs TDK vs RPA](https://rocketcea.readthedocs.io/en/latest/cea_tdk_rpa.html)
- Qualitative: "Low pressures: Kinetics; Low Thrusts: Boundary layer; Small Volume: Low combustion efficiency." Modern engines reach 95–98% Isp efficiency overall — [Karabeyoglu, Stanford AA284a Lecture 12](https://web.stanford.edu/~cantwell/AA284A_Course_Material/Karabeyoglu%20AA%20284A%20Lectures/AA284a_Lecture12.pdf)

### Inferences
- For rockmoc's default (where no ODK run is available), a fracKin regression on Rthrt, Pc, Tc and ε (the RocketIsp approach) is the cheapest upgrade. RocketIsp is open source, so its trained fracKin models or tables could be reused if the licence fits; check the licence.
- Background figures [RECALLED, UNVERIFIED]: large booster engines (F-1, RS-25, RD-180) have kinetic losses of about 0.1–0.5% Isp. Small storable (N2O4/MMH) thrusters of roughly 0.1–1 kN at low Pc lose about 1–3%. LOX/CH4 falls between RP-1 and LH2. Confirm these against SP-8120 and RocketIsp's TDK-derived plots.

### Gaps
- No tabulated η_kin versus Pc and thrust from JANNAF, SP-8120 or Huzel & Huang was retrieved. SP-8120 is at https://ntrs.nasa.gov/search.jsp?R=19770009165 (cited by RocketIsp), JPL TM 33-548 (simplified correlation procedures) is at https://ntrs.nasa.gov/citations/19730012958, and NASA CR-179025 (1986) is at https://ntrs.nasa.gov/search.jsp?R=19870009172. All three are worth fetching manually.

## 4. Two-phase (particle) losses: Hoglund, Kliegel, SPP, Hermsen D43

### Takeaway
Two-phase loss comes from velocity and thermal lag of condensed Al2O3 droplets. Larger droplets, higher particle loading, smaller throats and lower Pc all increase it. For aluminized solids the loss is a few percent of Isp. The JANNAF Solid Performance Program (SPP) computes it either empirically or with a 2-D two-phase MOC, using a Hermsen D43 particle size that depends on throat diameter. I could not retrieve the exact correlation coefficients.

### Cited Findings
- Mechanism: thermal lag ("the heat of the condensed-phase is partly ejected out of the nozzle before transferring this heat") and velocity lag ("particles exit the nozzle at a lower velocity than the gases"). The no-lag modified-γ equations "represent an upper limit on performance." Hoglund, "Recent Advances in Gas-Particle Nozzle Flows" (ARS J., May 1962), is cited for the isentropic-exponent treatment. Example: a 100 lbf motor "suffers a 6% loss in Characteristic Velocity if the average particle size is 1.5 micron" — [Nakka, Two-Phase Flow](https://www.nakka-rocketry.net/th_2phf.html)
- "Two phase flow loss is a result of the thermal and velocity lag between the expanding gas and particles." RocketIsp has no model and accepts a user-supplied η_TP — [RocketIsp](https://rocketisp.readthedocs.io/en/latest/nozzle_eff.html)
- SPP predicts performance through a set of independent efficiency models and includes two-phase effects. A user paper reports SPP empirical methods for kinetic, throat-erosion and submergence losses agreeing within 1% on four nozzles. This is a secondary source — [NWPU paper record](https://pure.nwpu.edu.cn/en/publications/calculation-of-nozzle-efficiency-for-solid-rocket-motor/); [SPP product description](https://archive.irixnet.org/siliconsurf/products/appsdirectory.dir/Applications/Computational_Fluid_Dynamics_CFD/ApplicationNumber6849.html)
- The 1990 improved-SPP paper treats D43 as "a unique function of motor throat diameter." As a rule of thumb, "1μ of droplet diameter … for every inch of nozzle throat diameter." An AFRPL finding, quoted secondhand, puts the Isp penalty at about 1% per 5% of added aluminum — [US patent 8051640 text](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8051640)
- At fixed inert-particle mass content, two-phase loss *decreases* with increasing particle size in a 2019 lab-scale-validated study. This contradicts the usual lag argument for aluminized motors; the cause is probably the regime (fine inert particles near equilibrium). Treat it as regime-specific — [BIT, Effects of particle size on two-phase flow loss](https://pure.bit.edu.cn/en/publications/effects-of-particle-size-on-two-phase-flow-loss-in-aluminized-sol)
- The primary report "Specific Impulse Losses in Solid Propellant Rockets" (AD-A004666) exists on DTIC but returned 403 to the fetch tool — [DTIC PDF](https://apps.dtic.mil/sti/tr/pdf/ADA004666.pdf)

### Inferences
- **1-D two-phase ODE [DERIVED, standard].** Add particle velocity u_p and temperature T_p, with particle loading as mass ratio φ = ṁ_p/ṁ_g:
```
du_p/dx = (u − u_p) / (τ_v u_p) · f_D        τ_v = ρ_p d² / (18 μ)        (Stokes; f_D = C_D Re_p/24 correction)
dT_p/dx = (T − T_p) / (τ_T u_p) · (Nu/2)     τ_T = ρ_p c_s d² / (12 k)   (Nu = 2 limit)
gas momentum: ρ u du/dx = −dp/dx − ρ φ u du_p/dx · (u_p-weighted; per unit gas flow)
gas energy:   ṁ_g (h + u²/2) + ṁ_p (c_s T_p + u_p²/2) = const
```
  Couple these to the gas 1-D equations, either frozen-γ or from section 1. Isp = (ṁ_g u_e + ṁ_p u_pe + p_e A_e)/(ṁ_g+ṁ_p)/g0, and η_TP = Isp / Isp_equilibrium-two-phase (no-lag; CEA with condensed Al2O3 gives that limit). This is Kliegel-type gas-particle nozzle flow. Kliegel's "constant fractional lag" approximation (u_p = K u, T_p − T_p∞ = L(T−…)) gives closed forms, but I did not retrieve the exact forms.
- **SPP empirical two-phase loss [RECALLED, UNVERIFIED — structure only].** I recall a form ΔIsp_2P (%) = C1·ξ^C2·D43^C3 / (Pc^0.15 · ε^0.08 · Dt^C4), with ξ = condensed-phase moles per 100 g, Pc in psia, Dt in inches, D43 in µm, and coefficient sets switching on Dt and D43 ranges. The exponents and C-values must be verified in Nickerson et al., SPP manual (AFRPL-TR-83-036) or Coats/Dunn papers. A targeted search did not confirm the 0.15 and 0.08 exponents.
- **Hermsen D43 [RECALLED, UNVERIFIED].** Hermsen (1981, J. Spacecraft Rockets 18(6), "Aluminum oxide particle size for solid rocket motor performance prediction") gives D43 [µm] = 3.6304·Dt^0.2932·[1 − exp(−0.0008163·ξ·Pc·τ_res)], with Dt in inches, Pc in psia and τ_res chamber residence time in s. Verify before use. The patent's "1 µm per inch of throat" rule is a rough consistency check.
- Typical loss [RECALLED, UNVERIFIED]: about 2–5% Isp for 16–20% Al composites. It is larger for small motors (Dt < 2–3 in) and smaller for large boosters (Dt > 20 in).

### Gaps
- The exact Hoglund correlation, the Kliegel closed-form loss, the SPP empirical constants and the Hermsen coefficients were all unverified. The original PDFs could not be fetched: AD-A004666 returned 403, and the AIAA/ARS papers are paywalled.

## 5. JANNAF performance-prediction methodology: combining divergence, BL, kinetic and two-phase losses

### Takeaway
Liquids follow the JANNAF/CPIA 246 approach: TDK, a 2-D MOC with finite-rate chemistry, captures divergence and kinetics together, including their interaction, and BL (TBL/BLM) is subtracted as a thrust decrement. In simplified methods the losses are multiplied as efficiencies. SPP for solids uses independent efficiency models, each evaluated as a loss relative to ideal. The BL term can be applied either as a factor or as a subtractive fraction, depending on how the correlation was built.

### Cited Findings
- Product form: Isp_del = Isp_ODE · η_ML · η_Vap · η_HL · η_FFC · η_Div · η_TP · η_Kin · η_BL — [RocketIsp efficiencies](https://rocketisp.readthedocs.io/en/latest/_sources/efficiencies.rst.txt)
- Alternate form: IspDel = IspODE · [ effML·effVap·effHL·effFFC·effDiv·effTP·effKin − (1 − effBL) ]. "Either form may be correct" depending on how the BL correlation was formulated. For preliminary design, take combined chamber losses (η_ERE) ≈ 0.98 — [RocketIsp efficiencies](https://rocketisp.readthedocs.io/en/latest/_sources/efficiencies.rst.txt)
- TDK and the BL module are coupled ("TDK/BLM") to produce a thrust-chamber Isp prediction — [Nickerson 1984 progress report PDF](https://ntrs.nasa.gov/api/citations/19860019159/downloads/19860019159.pdf)
- SPP computes "deviations from ideal performance with a set of independent efficiency models", covering kinetics, two-phase, divergence, BL, submergence and erosion — [SPP description / NWPU secondary](https://pure.nwpu.edu.cn/en/publications/calculation-of-nozzle-efficiency-for-solid-rocket-motor/)

### Inferences
- For rockmoc:
  - η_div comes from the existing MOC or a 2-D frozen/equilibrium integration.
  - η_kin comes from ODK (section 1) or Bray (section 2).
  - η_TP comes from the 1-D two-phase ODE or an SPP-type correlation.
  - η_BL comes from the existing BL estimate.
  - Combine as Isp = Isp_ODE·η_div·η_kin·η_TP − ΔIsp_BL, which is the subtractive-BL form above.
- Interaction terms. Divergence and kinetics do not separate exactly, because kinetic freezing depends on the streamline expansion rate, which differs on-axis and near the wall. Only a 2-D kinetic MOC (TDK) captures this. A rockmoc approximation: run ODK along the axis and along the wall streamline, then mass-average. Two-phase and divergence also interact, since particles do not follow gas streamlines and leave a particle-free zone near the wall. That interaction needs a 2-D two-phase MOC, as in SPP.

### Gaps
- The CPIA 246 text, with its exact JANNAF definitions of ΔIsp terms and how ODK, TDK and TBL are combined, was not retrieved. It is a controlled distribution document.
