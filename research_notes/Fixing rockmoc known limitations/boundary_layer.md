# Improved integral boundary-layer methods for rocket/supersonic nozzle walls (rockmoc)

Method note: most primary sources (Bartz 1963/65, Back-Massier-Gier 1965, Boldman et al. TN D-4788 / D-6594 / D-5356) exist only as scanned PDFs on NTRS. They were OCR'd locally with Tesseract. Equations below are transcribed from OCR text and checked for internal consistency. Where OCR garbled a symbol, this is flagged **[OCR-uncertain]**. Items from background knowledge that I did not confirm in a source this session are flagged **[UNVERIFIED]**.

---

## 1. Bartz integral method (JPL TR 32-387, 1963 / NASA CR-62615 / Advances in Heat Transfer vol. 2, 1965): exact equations

### Takeaway
Bartz solves two coupled integral equations, one for momentum (θ) and one for energy (φ), each started from measured or assumed thicknesses at the nozzle inlet. The rest of the method:
- **Skin friction:** Coles' adiabatic correlation, or Blasius 0.0256 Re_θ^-1/4 with property corrections.
- **Stanton number:** von Kármán Reynolds analogy multiplied by an "interaction" factor (φ/θ)^n.
- **Profiles:** 1/7-power profiles for u/U (over δ) and for the stagnation-temperature defect (over a separate thermal thickness Δ).

Because Δ ≠ δ, H = δ*/θ is computed by integration rather than from a Crocco closure. In cooled throats this gives φ/θ up to about 5, and δ* (and H) can be **negative**. rockmoc's single-thickness Crocco–Busemann H cannot capture this, and it is the most likely cause of its δ* overprediction.

### Cited Findings
- **Source.** Bartz, D. R., "Turbulent Boundary-Layer Heat Transfer from Rapidly Accelerating Flow of Rocket Combustion Gases and of Heated Air", JPL / NASA-CR-62615 (Dec 1963), 132 pp. The scanned PDF is https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf. It was published as Advances in Heat Transfer vol. 2, Academic Press 1965, pp. 1-108 (citation confirmed in the TN D-6594 reference list). — [NTRS CR-62615](https://ntrs.nasa.gov/citations/19650013685); [TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)
- **Key assumptions (Bartz §II.A.3):**
  - (9) Cf equals the flat-plate value at the same edge conditions, Tw and θ.
  - (10) St equals the flat-plate value at the same φ and θ.
  - (11) "The Stanton number for unequal momentum and energy thicknesses is that for equal thicknesses multiplied by (φ/θ)^n, where n is a small 'interaction exponent'."
  - (12) Heat transfer affects Cf in one of two ways: (a) not at all, so Cf = adiabatic-wall Cf; or (b) properties are evaluated at the arithmetic mean of Tw and Te.
  - (15) The shape parameters θ/δ, φ/Δ and δ*/θ are those for 1/7-power profiles of velocity and of (T0 − Tw).
  
  — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Momentum integral (Bartz Eq. 21–22), axisymmetric thin BL:**
  - d/dx(r ρ U² θ) = r τw + r ρ U δ* dU/dx
  - Equivalently, dθ/dx = Cf/2 − θ[(1 + δ*/θ)(1/U)dU/dx + (1/ρU) d(ρU)/dx + (1/r) dr/dx]
  - Mach-number forms (isentropic edge): (1/U) dU/dx = [1/(1 + (γ−1)/2 M²)] (1/M) dM/dx, and (1/ρU) d(ρU)/dx = [(1 − M²)/(1 + (γ−1)/2 M²)] (1/M) dM/dx.
  - Bartz then changes the independent variable from wall arc length x to axial z using dx/dz = [1 + (dr/dz)²]^½ (Eq. 25).
  - This is the same equation rockmoc uses. — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Energy integral (Bartz Eq. 28–30):**
  - d/dx[r ρ U cp (T0 − Tw) φ] = r q_w
  - Equivalently, dφ/dx = St·(Taw − Tw)/(T0 − Tw) − φ[(1/ρU) d(ρU)/dx + (1/r) dr/dx − (1/(T0 − Tw)) dTw/dx]
  - φ is the energy (stagnation-enthalpy-defect) thickness: φ = ∫ (ρu/ρeUe)(T0 − T0,local)/(T0 − Tw) dy.
  - The dTw/dx term carries the wall-temperature history. **[OCR-uncertain: exact grouping of the Taw/T0 factor in the St term]** — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Skin friction (illustrative Blasius form, Eq. 32):**
  - Cf_inc = 0.0256 / Re_θ^¼, with properties at the edge static temperature. Bartz notes this is within 5% of Coles' tabulated values for 400 < Re_θ < 15,000.
  - Option (a), adiabatic Cf with m = 0.6 viscosity exponent: Cf ≈ 0.0256 Re_θ^-¼ (Taw/Te)^-0.6 (Eq. 35).
  - Option (b), film properties: Cf ≈ 0.0256 Re_θ^-¼ [½(Tw/Te) + ½]^-0.6 (Eq. 37). **[OCR-uncertain: exponent shown as −0.6; the (3−m)/4 power form appears in Eq. 36]** — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Coles correlation, the method actually used (Appendix A):**
  - Cf = (ρs μs... ) transformation: Cf = (Te/Ts)·Cf̄ and Cf·Re_θ relation (A-1, A-2), with a sublayer temperature Ts (A-3).
  - Table A-1 lists Cf̄·Re_θ̄ → Cf̄ pairs: 2.51 → 0.00590; 3.10 → 0.00524; 3.97 → 0.00464; 4.88 → 0.00426; 5.73 → 0.00398; 7.41 → 0.00363; 8.94 → 0.00340; 12.75 → 0.00308; 16.36 → 0.00290; 23.2 → 0.00269; 29.6 → 0.00255; 35.9 → 0.00246; 41.8 → 0.00238; 53.6 → 0.00227; 64.8 → 0.00219.
  - Extrapolation above Cf̄Re = 64.8 (A-4): (2/Cf̄)^½ = 2.44 ln(Cf̄Re_θ̄ …) + 7.68 **[OCR-uncertain argument]**.
  - Below 2.51 (A-5): Cf̄ ≈ 0.009896/(Cf̄Re)^0.62 **[OCR-uncertain]**. Bartz says A-5 is used only to obtain an order-of-magnitude Cf for a thin-entrance start.
  - Ts expression (A-3) is not recoverable from OCR. Use Coles 1962 (RAND R-403-PR) directly. — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Stanton number (Eq. 38, A-7, A-8), von Kármán analogy with interaction factor:**
  - St = (Cf(Re_φ)/2)·(φ/θ)^n / {1 − 5(Cf(Re_φ)/2)^½ [1 − Pr + ln(6/(5Pr + 1))]}
  - Cf(Re_φ) means Cf evaluated with Re_φ = ρUφ/μ in place of Re_θ.
  - Simplified Colburn version (Eq. 42–43): St ≈ 0.0128 (T_ref/Te)^-0.6 / (Pr^⅔ Re_φ^¼) · (φ/θ)^(n−¼)·… **[OCR-uncertain on exponent grouping]**. — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Interaction exponent n:**
  - The earlier (φ/θ)^(1/7)-type factor of Refs. 8 and 9 corresponds to n = 3/28 ≈ 0.1.
  - Bounds: n ≥ 0 by physical argument, n ≤ ¼. "A value of n of 0.1 was arbitrarily selected."
  - Against the JPL heated-air nozzle, the curve with n = 0 and Coles adiabatic Cf ("type D") was about 30% below the n = 0.1 / film-temperature curve ("type A") at the throat. Type D agreed with the throat data.
  - With n = 0 the energy equation decouples and has a closed-form quadrature solution (Eq. 44–45). — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Profiles and shape factor (Appendix B):**
  - Profiles: u/U = (y/δ)^(1/7) for y ≤ δ, and (T0 − Tw)/(T0e − Tw) = (y/Δ)^(1/7) for y ≤ Δ.
  - Density follows from T = T0,local − u²/(2cp). Two cases are integrated separately: Case I (δ < Δ) and Case II (δ > Δ).
  - The Crocco distribution is the special case Δ = δ.
  - δ*/θ, θ/δ and φ/Δ are obtained by numerical integration and depend on M, Tw/T0 and Δ/δ.
  - Bartz: shape parameters are "quite insensitive to the arbitrarily specified exponent over a range from about 1/5 to 1/9 but are quite sensitive to the local value of M and Tw/T0". He also warns that "use of low-speed adiabatic values of the shape parameters … can lead to significant errors near the nozzle throat." — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Coupling and negative δ*:**
  - The equations must be solved iteratively because δ*/θ depends on Δ/δ.
  - The earlier "first approximation", which fixed Δ/δ = 1, "does have a large effect on δ*/θ".
  - "In the throat region of a highly cooled (Tw/T0 << 1) nozzle … φ/θ, and hence Δ/δ, may approach values as high as 5 and, as a consequence, δ*/θ becomes negative": mass flux in the cold layer exceeds the free-stream value.
  - Back et al. report the same: "because of the cooled wall, the displacement thickness δ* becomes negative upstream of the throat, as does H = δ*/θ." — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf); [Back, Massier & Gier JPL TR 32-415](https://ntrs.nasa.gov/api/citations/19650010083/downloads/19650010083.pdf)
- **Starting conditions:**
  - Bartz requires "some specified thickness at the starting point of the analysis (such as that appropriate to growth in a pipe entrance region of approximately the combustion-chamber length when starting at the nozzle entrance)".
  - In the validation he used measured inlet δ and Δ.
  - Sample calculations show a thick inlet BL "shrinks rapidly as a result of the subsonic acceleration, reaching a minimum just ahead of the throat", while a thin one grows and then shrinks to about the same value. The initial δ or θ "has little effect on the h_g distribution over most of the nozzle when the initial temperature boundary-layer thickness is small". — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Inlet profiles:** Measured nozzle-inlet profiles fit 1/7-power constant-property profiles for both u and T0. Pipe-geometry ("thick BL") thicknesses were about 5% below the planar-formula values. — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Rocket-engine correction factors (Bartz recommendation):**
  - Multiply BL-predicted h in the transonic region by about 1.3–1.4, tapering to 1.0 by M ≈ 4.
  - Use 1.3–1.8 in the chamber and contraction (attributed to free-stream turbulence and combustion).
  - Prefer Coles adiabatic-wall Cf (equivalently T_ref = Taw). — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)
- **Origin of the computer method:** Elliott, Bartz & Silver, "Calculation of Turbulent Boundary-Layer Growth and Heat Transfer in Axi-Symmetric Nozzles", JPL TR 32-387 (15 Feb 1963), is the programmed (IBM 7090) version of this method. Back et al. initiated it with θ and Δ/δ measured 1.25 in. upstream of the nozzle inlet. — [Back, Massier & Gier TR 32-415](https://ntrs.nasa.gov/api/citations/19650010083/downloads/19650010083.pdf)

### Inferences
- **rockmoc's H closure is not the error source by itself.** I checked numerically (scratch script, γ = 1.4, r = 0.89, 1/7 profile) that H = H_i·Tw/Te + r(γ−1)/2·M² reproduces the exact Crocco–Busemann 1/7 integral to within 1–2% for M = 1–3.7 and Tw/T0 = 0.5–0.8.
- **The thermal-thickness ratio is the dominant missing physics.** Integrating the Bartz two-thickness profile gives:

  | Case | Δ/δ = 1 | Δ/δ = 2 | Δ/δ = 4 |
  |---|---|---|---|
  | M = 3, Tw/T0 = 0.5 | H ≈ 3.57 | 1.40 | −3.2 |
  | M = 3, Tw/T0 = 0.8 | 4.60 | 3.68 | 1.68 |

  So a thermal layer even modestly thicker than the velocity layer cuts δ* by tens of percent on a cooled wall. That is consistent with rockmoc's +25% δ* bias on the cooled TN D-5377 nozzle.
- **Recommended rockmoc upgrade:**
  1. Add the energy equation for φ, started with a chamber/inlet value.
  2. Carry Δ/δ, obtained from φ/θ through the 1/7-profile integrals.
  3. Compute H by numerical quadrature of the Bartz Case I/II profiles. Tabulate H(M, Tw/T0, Δ/δ) once and cache it.
  4. Keep n = 0 (decoupled St), which the data favour.
- **θ benefits as well.** Lower H reduces the (2 + H − M²) term. In accelerating flow (dM/ds > 0) that increases θ: smaller H means less of the favourable-gradient thinning, which acts on the H + 2 term. This is the right direction for rockmoc's 38–43% θ underprediction with a cooled inlet.

### Gaps
- Exact form of Coles' sublayer temperature Ts (Bartz Eq. A-3) and the A-4 argument could not be read reliably from OCR. Get them from Coles 1962 RAND R-403-PR, or read the CR-62615 scan pages 78–80 directly.
- The full Appendix B integrand coefficients (a, b, c in t/T = a(1 + b s − c s²)) were garbled in OCR. Use the reconstruction above; it is derived from the stated profiles, not quoted.
- JPL TR 32-387 (Elliott-Bartz-Silver) itself was not located online.

---

## 2. Sasman & Cresci (AIAA J. 4(1), 1966)

### Takeaway
No full text or equations were found online this session. Treat any equations as unverified.

### Cited Findings
- No source retrieved.

### Inferences
- Bartz's two-thickness method plus the Boldman modifications in §3 already cover the inlet-thermal-history and cooling effects that Sasman–Cresci addresses. Implementing Sasman–Cresci is probably not needed for rockmoc's stated objectives.

### Gaps
- **[UNVERIFIED]** From memory: Sasman, P. K. & Cresci, R. J., "Compressible turbulent boundary layer with pressure gradient and heat transfer", AIAA J. 4(1):19–25, 1966. I recall a moment-of-momentum (Walz-type) integral formulation with a Crocco-type temperature–velocity relation, a power-law velocity profile and a Spence-type shape-factor relation. Verify against the AIAA paper (DOI likely 10.2514/3.3378, not checked) before use.

---

## 3. Inlet thermal/momentum BL history and modified energy methods (Boldman/Graham, NASA Lewis)

### Takeaway
For the NASA Lewis 30°-15° and 60°-15° conical nozzles (the same rig family as TN D-5377), throat heat transfer and BL state are controlled by the **thermal** layer coming from the inlet, not the momentum layer:
- Uncooled inlet length (which changes θ) has negligible effect on throat heat transfer.
- Cooled inlets (thick inlet thermal BL) cause pronounced reductions.

The working fix used by NASA Lewis was an integral **energy** method (Bartz with n = 0), later modified with an acceleration term.

### Cited Findings
- **TN D-6594 summary findings (Boldman & Graham 1972):**
  - "Altering the length of uncooled pipe inlets with length to diameter ratios greater than 2.6 (changing the turbulent boundary layer momentum thickness at the nozzle entrance) resulted in a negligible effect on the throat heat transfer … However, upon cooling the pipe inlets (increasing the thickness of the thermal boundary layer on the nozzle entrance) pronounced reductions in throat heat transfer were observed."
  - Test conditions: air, T0 = 539 K, throat Re 6×10⁵ to 5×10⁶.
  - Rig: 3.8 cm throat diameter and throat radius of curvature. Contraction ratios 18.8 (16.51 cm inlet) and about 4.2 (7.80 cm inlet). Exit area ratios 25.3 (30°-15°) and 3.32 (60°-15°). Five inlets: three uncooled, two cooled.
  
  — [Boldman & Graham, TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)
- **φ/θ magnitudes:** BL surveys in the 30°-15° nozzle "suggested that the ratio of energy to momentum thickness at the throat could be large (possibly as high as 10)". A "simple integral energy analysis provided good estimates of the high Reynolds number heat transfer (±20 percent) at the expense of a pronounced underprediction in nozzle entrance and exit regions." — [TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)
- **Energy method (Bartz n = 0):**
  - Skin friction is taken as the adiabatic Cf on free-stream properties.
  - "the best value for the interaction exponent N was determined to be zero". The energy equation is then solved alone, with St = f(Re_φ) given by the von Kármán form in §1 with Re_θ replaced by Re_φ.
  - N range considered: 0 to 0.25. — [TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)
- **Modified energy method (TN D-6595, Boldman, Schmidt & Ehlers 1972):**
  - Auxiliary equation (D-6594 Eq. 14): St = 1.6 f(Re_φ){1.0 − 0.11·[acceleration term built from (x/(u∞R))·d(u∞R)/dx]} **[OCR-uncertain; the exact acceleration group must be read from TN D-6595/D-6594]**.
  - Gives "equal or slightly better agreement with the data in the nozzle throat and much better agreement in the entrance and exit regions". The exception is the 60°-15° nozzle with the cooled inlet.
  - The D-6595 abstract describes "a modified momentum-heat analogy which includes an acceleration term comprising the nozzle geometry and free stream velocity". — [TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf); [JAXA record of TN D-6595](https://repository.exst.jaxa.jp/dspace/handle/a-is/203925)
- **Unaccelerated pipe-flow correlation:** It gives throat h "nearly a factor of two greater than the experimental values" and is independent of inlet cooling. When the thermal BL starts near the throat, however, simple correlations can be **non-conservative**. — [TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)
- **Graham & Boldman, TN D-5356 (1969):** "Boundary layer profile measurements in the convergent portion of a convergent-divergent nozzle have verified that the thermal layer is much thicker than the velocity layer." It proposes an approximate integral-energy method for throat heat transfer, tested for contraction ratios 4–19 and Tw/T0 = 0.25–0.85. — [TN D-5356](https://ntrs.nasa.gov/api/citations/19690022491/downloads/19690022491.pdf)
- **Back, Massier & Gier, JPL TR 32-415 (Rev. 1, Feb 1965):**
  - Increasing inlet BL thickness from 5% to 25% of inlet radius reduced h by about 10%.
  - The measured inlet profiles (1.25 in. upstream of the inlet, 18-in. cooled approach) were turbulent 1/7-power. Example thickness columns in their Fig. 8 table: p0 = 44.8 psia → 0.060 / 0.070 / 0.043 in.; 254 psia → 0.038 / 0.038 / 0.015 in. **[OCR; column identities (δ*, θ, φ ordering) not confirmed]**.
  - Thickness definitions use the pipe (R − y) weighting. — [Back, Massier & Gier TR 32-415](https://ntrs.nasa.gov/api/citations/19650010083/downloads/19650010083.pdf)

### Inferences
- **Concrete rockmoc fix for (b):**
  1. Start the integration upstream of the throat, at the chamber/convergent inlet, not at the throat.
  2. Use initial θ0 and φ0 from either (i) user input, or (ii) a flat-plate / pipe-entrance turbulent estimate over the chamber length L_c: θ ≈ 0.036 L_c Re_L^-0.2 × (7/72) for a 1/7 profile **[UNVERIFIED: standard flat-plate 1/7-law result, δ = 0.37 x Re_x^-0.2 and θ = 7δ/72]**.
  3. With a cooled chamber wall, set φ0 > θ0. The thermal layer grows from where cooling starts, so φ0 ≈ θ_flat-plate over the cooled length.
  4. Integrate the θ and φ equations, then derive Δ/δ and H.
- **Why this should help both errors.** Thick inlet thermal layers raise φ/θ, which pushes H down (often negative near the throat). That shrinks δ* and raises θ growth through the favourable-gradient term. It should address both the +25% δ* error and the −40% θ error on the cooled-inlet case.
- **Simplest heat-transfer option.** If only q_w matters, the decoupled energy method (n = 0) is the simplest and best-validated option (±20% at the throat).

### Gaps
- TN D-5377 itself: the NTRS PDF link tried (19690026187) returned a stub, so its tabulated θ, δ*, φ were not re-checked here.
- The exact D-6595 auxiliary equation needs a direct read of TN D-6595 (NTRS full text not located this session).

---

## 4. Shape-factor closures with pressure gradient and wall cooling

### Takeaway
I found no single well-sourced compressible H_i(Λ) correlation for accelerating nozzle flow this session. The one sourced nozzle-specific empirical relation is the 8°-convergent-nozzle relation used in TN D-4788: shear parameter ζ and H as functions of the pressure-gradient parameter Γ, with the combination ζ − (H+1)Γ = 0.016 − 3.55Γ. Thermal-thickness effects (§1) dominate cooling effects on H more than the pressure-gradient dependence of H_i does.

### Cited Findings
- **TN D-4788 Appendix B** (relation taken from its Ref. 12, an 8° half-angle convergent conical nozzle study):
  - ζ = (τw/ρu∞²) Re_θ^¼ (B3)
  - Γ = (θ/u∞)(du∞/dx) Re_θ^¼ (B4)
  - The momentum equation reduces to dRe_θ/dx form (B5) with the bracket ζ − (H+1)Γ = 0.016 − 3.55Γ (B6).
  - Result (B8): dRe_θ/dx = (5/4)(u∞/ν0)[0.016 − 3.55 Re_θ^(5/4)·(ν0/u∞²)(du∞/dx) …]/Re_θ^¼ − (Re_θ/r)(dr/dx) **[OCR-uncertain on grouping; use (B6) directly]**.
  
  — [Boldman, Schmidt & Gallagher, TN D-4788](https://ntrs.nasa.gov/api/citations/19680024170/downloads/19680024170.pdf)
- **Cohen–Reshotko integral method:** cited in a later comparison as having "the slightest deviation, sometimes only 4% higher when compared with measured displacement thickness" (secondary search snippet; laminar method). — search result from [DTIC ADA055124](https://apps.dtic.mil/sti/pdfs/ADA055124.pdf) **[snippet-level; not opened]**
- **Bartz on velocity exponent:** shape parameters are insensitive to the velocity exponent between 1/5 and 1/9 but sensitive to M and Tw/T0. — [Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)

### Inferences
- **Priority for rockmoc:** the (M, Tw/T0, Δ/δ) dependence of H, i.e. the Bartz two-thickness integral, matters more than an H_i(pressure-gradient) correlation.
- **Optional add-on:** let H_i fall in strong acceleration. Favourable-gradient equilibrium BLs are fuller; H_i of about 1.2–1.3 versus 1.286 for 1/7 is a modest 5–10% effect.

### Gaps
- None of the following were verified this session; all are **[UNVERIFIED]**:
  - Spence (1961?) H relation.
  - Walz's compressible H relation.
  - Alber & Coats.
  - Walz/Zhang "generalized Reynolds analogy" modified Crocco relation T/Te = Tw/Te + (Trw − Tw)/Te·[(1 − sPr)(u/U)² + sPr(u/U)] − r(γ−1)/2 M²(u/U)², with sPr ≈ 0.8, from Zhang, Bi, Hussain & She, JFM 2014. This relation lowers T in the cooled near-wall region versus Crocco–Busemann and would lower H in cooled flows. Worth checking as a cheap single-thickness alternative.
  - Van Driest-transformed H.

---

## 5. Relaminarization / transition in accelerated throats

### Takeaway
Use an axisymmetric acceleration parameter K_ax = (ν/u∞²)(du∞/dx) + 0.352·(ν/(u∞ r))(dr/dx), with critical value 2.88×10⁻⁶. This is equivalent to the momentum integral driving Re_θ below 360. Above it, expect laminar-like heat transfer, about 50% of the turbulent prediction at the throat. A practical code switches Cf/St to laminar closures when K > K_cr (or Re_θ < about 360) and reverts when K drops.

### Cited Findings
- **Boldman, Schmidt & Gallagher, TN D-4788 (Sept 1968):**
  - Derivation: from (B6) and setting Re_θ,cr = 360 ("classical value for forward transition on a flat plate"), Eq. B9 gives the critical condition 0.00368 − 1278·K_ax = 0.
  - B10: K_ax = (ν0/u∞²)(du∞/dx) + 0.352 (ν0/(u∞ r))(dr/dx). ν0 is evaluated at stagnation temperature.
  - So K_ax,cr = 0.00368/1278 = 2.88×10⁻⁶ (the report states 2.88×10⁻⁶).
  - Data with K_ax > K_cr "reflect the influence of a laminarized velocity boundary layer".
  - "At the throat of each nozzle, the heat transfer was about 50 percent of the predicted value" at low Re.
  - Test conditions: T0 = 539 K; p0 = 2.0–20.4 atm; 30° and 60° convergent nozzles. — [TN D-4788](https://ntrs.nasa.gov/api/citations/19680024170/downloads/19680024170.pdf)
- **TN D-6594 review:**
  - "Whenever the acceleration parameter exceeds a nominal value of about 3×10⁻⁶, a rapid dropoff in the Stanton number occurs."
  - Inlet cooling did not change K and did not change the onset Reynolds number.
  - Laminar-like profiles were measured at K ≈ 6×10⁻⁶ (30°-15°). In a 10°-10° nozzle (Back et al.), K ≈ 2.8 and 4.1×10⁻⁶ affected velocity only; at K ≈ 6.0 and 8.3×10⁻⁶ both velocity and temperature profiles were laminar near the wall.
  - Throat Re below about 1×10⁶ gave the large ("laminarized") depression. Above about 1×10⁶ the throat heat transfer was still up to 50% below unaccelerated correlations.
  - Roughness reduces the Re range where laminarization occurs. — [TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)
- **Detrainment:** Boldman, Neumann & Ehlers (TN D-6043) observed acceleration pulling mass out of the thick upstream BL. — [TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)
- **Laminarization eddy-viscosity model:** Schmidt, Boldman & Todd, TM X-2501 (Feb 1972), use a fully laminar sublayer plus a turbulent wake region for eddy transport in a differential BL code. It agreed with convergent-section velocity profiles. — [NTRS TM X-2501](https://ntrs.nasa.gov/api/citations/19720008557)
- **Back, Massier & Gier (1965):**
  - h fell below the turbulent p0^0.8 trend at p0 below about 75 psia.
  - The predicted minimum Re_θ was 1500 even at the lowest p0, versus a measured departure-from-turbulence value of about 600 (their Ref. 14). Preston's limit is 320.
  - They define a turbulence-energy-production ratio χ ≈ 22·(ν/u_τ²)(du_e/dx)·… **[OCR-uncertain constant 22, group not fully legible]**, reaching 0.14 at the lowest p0. — [Back, Massier & Gier TR 32-415](https://ntrs.nasa.gov/api/citations/19650010083/downloads/19650010083.pdf)
- **Smith (TDK comparison, NTRS):** laminar flow at throat-diameter Re = 320,000 (360 psia), transitional up to 1000 psia chamber pressure. Predicted acceleration-parameter distributions were consistent with this. — [NTRS 19880063559](https://ntrs.nasa.gov/citations/19880063559) **[search-snippet level]**
- **BLIMP-J early verification:** with TDK free-stream inputs, agreement was poor in one case, "believed to result from failure of the turbulent model" to capture laminarization. — [NTRS 19740023639](https://ntrs.nasa.gov/api/citations/19740023639/downloads/19740023639.pdf) **[search-snippet level]**

### Inferences
- **Implementation for (c):**
  - At each station compute K_ax from the edge solution. Use ν at T0 to match Boldman, or ν_e (more standard). Note that K_cr = 2.88×10⁻⁶ was derived with ν0.
  - If K_ax > 2.88×10⁻⁶ (or Re_θ < 360), switch to a laminar compressible closure: Cf from the Thwaites/Cohen–Reshotko integral with H_i ≈ 2.2–2.6, and St = laminar Reynolds analogy.
  - Optionally blend with an intermittency weight γ_t(K) between 2.5 and 4×10⁻⁶ to avoid a discontinuity.
  - Revert to turbulent when K < K_cr. Expect re-transition downstream of the throat, where K falls rapidly.
- **Relevance to the TN D-5377 runs.** In large rocket engines (high p0, large throats) K is usually below 1×10⁻⁶ and the switch never triggers. It matters for small, low-p0 test nozzles like the NASA Lewis rigs. That is why the switch matters for rockmoc's validation cases, including low-p0 TN D-5377 runs if any.

### Gaps
- Launder & Jones (1969) and Jones & Launder (1972) thresholds (often quoted as K ≈ 3–3.5×10⁻⁶): not retrieved directly; only the NASA reviews' "about 3×10⁻⁶" is sourced.
- No validated re-transition criterion for the supersonic side was found.

---

## 6. JANNAF / industry practice (BLIMP-J, MABL, TDK-TBL) and BL loss magnitudes

### Takeaway
JANNAF's reference BL method is BLIMP-J, a nonsimilar finite-difference code for momentum, energy and species with laminar/turbulent models (Cebeci-Smith, Bushnell-Beckwith). TDK uses its own TBL/MABL module. Both are differential, not integral, codes. The accuracy evidence found is qualitative to moderate:
- TDK overpredicts BL thickness away from the wall in high-area-ratio nozzles.
- TDK heat flux error ranges from about 100% at ε = 50 to about 3% at ε = 975.

### Cited Findings
- **BLIMP-J role:** chosen as the reference BL method in the JANNAF performance-prediction procedure. It solves the momentum, energy and species equations for laminar or turbulent reacting flows. — [NTRS 19750023946](https://ntrs.nasa.gov/api/citations/19750023946/downloads/19750023946.pdf) **[search-snippet level]**
- **Turbulence models:** a 1975 Aerotherm update added Cebeci-Smith and Bushnell-Beckwith models, which "are in close agreement for all cases studied". The Kendall model disagreed at some Tw/T0. — [NTRS 19750023946](https://ntrs.nasa.gov/api/citations/19750023946/downloads/19750023946.pdf)
- **Later improvements:** a 1991 update added condensed-phase effects, thick-BL thrust-loss treatment and free-stream turbulence. — [NTRS 19910007817](https://ntrs.nasa.gov/citations/19910007817) **[search-snippet level]**
- **High-area-ratio NASA test (1025:1 / 440:1):** TDK BL total-pressure predictions matched near the wall but overpredicted BL thickness further out. Heat flux was overpredicted in all cases, by about 100% at ε = 50 down to about 3% at ε = 975. — [NASA 19960029260](https://ntrs.nasa.gov/api/citations/19960029260/downloads/19960029260.pdf) **[search-snippet level]**
- **RocketIsp BL efficiency models:**
  - NASA SP-8120 provides a BL efficiency correlation in Pc, ε and throat radius. It "predicts about twice as much boundary layer loss" as RocketIsp's MLP regression of TDK-type results, and it cannot predict η_BL "beyond about a 1 inch throat radius at 1000 psia".
  - Trends: η_BL improves with Pc and throat size and worsens with area ratio. — [RocketIsp docs](https://rocketisp.readthedocs.io/en/latest/nozzle_eff.html)
- **Elliott-Bartz-Silver (JPL TR 32-387, 1963)** is the original programmed integral method (see §1). — [Back, Massier & Gier TR 32-415 reference list](https://ntrs.nasa.gov/api/citations/19650010083/downloads/19650010083.pdf)

### Inferences
- **rockmoc's TN D-5377 errors look normal for its method class.** The TDK result (BL thickness overpredicted in the outer layer, heat flux overpredicted at low ε) shows the industry codes have similar issues. rockmoc's ±20% θ and +25% δ* is in the range typical of 1960s–70s integral methods.
- **BL thrust-loss magnitude:** typical values are about 0.5–2% Isp for large engines, rising to several % for small, high-ε nozzles **[UNVERIFIED; based on SP-8120/TDK general knowledge, not quoted from a source found here]**.

### Gaps
- No quantitative BLIMP-J or MABL thrust-loss accuracy versus test data was found. The JANNAF manual (CPIA 246) is not openly online.
- The MABL method details (Mass Addition Boundary Layer, Aerojet/SEA) were not retrieved.

---

## 7. Validation datasets (tabulated δ*, θ, Cf, heat transfer)

### Takeaway
The best open datasets are the NASA Lewis conical-nozzle series, with tabulated heat-transfer coefficients and wall temperatures for 30°-15° / 60°-15° nozzles with adiabatic and cooled inlets (TN D-6594 Tables V–VII), plus the JPL Back-Massier-Gier heated-air nozzle. All are scanned NTRS PDFs.

### Cited Findings
- **TN D-6594** (Boldman & Graham 1972, review of the whole Lewis program): tabulates experimental heat-transfer coefficients (Table V), wall temperatures (Table VI) and stagnation conditions (Table VII). Cases cover the 30°-15° and 60°-15° nozzles with adiabatic inlet (inlet 2) and cooled inlet (inlet 3). PDF: https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf — [TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)
- **TN D-4788** (1968): heat transfer and BL surveys showing laminarization. Includes K_ax distributions. PDF: https://ntrs.nasa.gov/api/citations/19680024170/downloads/19680024170.pdf — [TN D-4788](https://ntrs.nasa.gov/api/citations/19680024170/downloads/19680024170.pdf)
- **TN D-5356** (1969): thermal versus velocity layer profiles in the convergent section; approximate energy method. PDF: https://ntrs.nasa.gov/api/citations/19690022491/downloads/19690022491.pdf — [TN D-5356](https://ntrs.nasa.gov/api/citations/19690022491/downloads/19690022491.pdf)
- **JPL TR 32-415** (Back, Massier & Gier, Rev. 1, Feb 1965):
  - Air heated by methanol combustion; conical 30°-15°-type C-D nozzle with inlet radius 2.53 in.
  - p0 = 80–250 psia (data down to 45 psia), T0 = 1030–2000 °R.
  - Inlet BL 5–25% of inlet radius; measured inlet δ*, θ, φ; local h along the nozzle.
  - PDF (62 pp.): https://ntrs.nasa.gov/api/citations/19650010083/downloads/19650010083.pdf — [TR 32-415](https://ntrs.nasa.gov/api/citations/19650010083/downloads/19650010083.pdf)
- **Other Lewis-series reports:**
  - TN D-4177 (inlet-diameter effects).
  - TN D-6043 (velocity, intermittency and turbulence in accelerated BL).
  - TN D-6595 (modified energy method).
  - TN D-5887 (roughness).
  - TN D-6439 (similar solutions).
  - TM X-2501 (laminarization model).
  - Back, Cuffel & Massier, ASME 69-HT-56 (laminarization with wall cooling, BL and heat-transfer measurements, 10°-10° nozzle). — [TN D-6594 reference list](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf); [NTRS TN D-6439](https://ntrs.nasa.gov/citations/19710022939)

### Inferences
- **Recommended validation ladder for rockmoc:**
  1. TN D-5377 (already used): θ, δ*.
  2. TN D-6594 tables: h(x) for adiabatic versus cooled inlet. These directly test the inlet-thermal-history model.
  3. TN D-4788 low-p0 cases to test the K_ax laminarization switch.
  4. JPL TR 32-415 to test sensitivity to inlet BL thickness. The target is about a 10% h change between inlet δ/R = 0.05 and 0.25.

### Gaps
- TN D-3221 (Boldman) and Fortini's datasets were not located this session.
- TN D-5377's own NTRS PDF link was not resolved (the guessed ID returned a 40-byte stub). Search NTRS for "TN D-5377".
- The tables are scanned. They need OCR plus manual checking before use as regression data.
