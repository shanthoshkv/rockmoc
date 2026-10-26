# Reference data for validating axisymmetric MLN and ideal-contour MOC (rockmoc)

Conventions used below: y* or r* = throat radius (not diameter) unless stated. "GP" = calorically perfect gas, "HT" = Zebbiche's high-temperature calorically imperfect model (Cp(T), γ(T)). All transcribed numbers are copied from the extracted PDF or text of the source. Values I derived myself are labelled **[derived]**.

## Q1. Argrow & Emanuel (1988) tables, and reproductions of axisymmetric MLN numbers

### Takeaway
I could not get the Argrow & Emanuel (1988) table. The paper is paywalled, and no open-access source I found reproduces its axisymmetric numbers. The best open, citable tabulation of axisymmetric sharp-corner MLNs is T. Zebbiche's Doctorat d'État thesis (Univ. Blida). It has a GP column (γ = 1.402) for M_S = 1.5 to 6, giving L/y*, corner angle θ*, max wall angle θ_max, exit radius, and the inflection point. These rows can be used directly as a regression table. The 2024 INCAS paper by Abada et al. reproduces the HT (T0 ≈ 2000 K) rows of that thesis, not the GP rows.

### Cited Findings
- **Argrow & Emanuel citation.** Argrow BM, Emanuel G (1988), "Comparison of minimum length nozzles", J Fluids Eng 110(3):283–288, doi:10.1115/1.3243546. A 2025 open-access paper cites it as a reference for MLN design, but the citing paper does not reproduce its numbers. — [Advances in Aerodynamics 2025, s42774-025-00244-1](https://link.springer.com/article/10.1186/s42774-025-00244-1)
- **Zebbiche thesis: what it covers.** Title: "Conception et dimensionnement des tuyères supersoniques de propulsion", T. Zebbiche, Thèse de Doctorat d'État, Université Saad Dahleb de Blida. It covers 2D and axisymmetric "tuyère à détente centrée" (MLN), GP and HT, solved by MOC with finite differences and a predictor-corrector scheme. The abstract states: "axisymmetric nozzle has an inflexion point which is not the case for the bidimensionnal nozzle." — [Blida repository PDF](https://di.univ-blida.dz/jspui/bitstream/123456789/6224/1/32-620-91-1.pdf)
- **Zebbiche thesis: definitions.** These come from §3, near the 2D wall-point construction, and are reused for the axisymmetric tables. The wall passes through the corner point A, where x_A = 0, y_A = y*, M_A = M*, θ_A = θ*. M* is "the Mach number at point A just after the expansion", i.e. the first point on the last descending characteristic AB. θ* is the flow angle at the corner (the expansion angle). θ_max is the largest wall angle. In the axisymmetric case θ_max occurs at the downstream inflection point (Tables 5.48/5.49), not at the corner. — [same thesis](https://di.univ-blida.dz/jspui/bitstream/123456789/6224/1/32-620-91-1.pdf)
- **Zebbiche thesis, Tables 5.42–5.47: axisymmetric MLN, GP rows (γ = 1.402).** Columns are, in order: L_d/y* (Kernel length; the thesis's Fig. 5.71 is titled "longueur de Kernel"), L/y* (nozzle length), θ* (°), θ_max (°), M*, y_S/y* (exit radius / throat radius), C_Masse (nondimensional nozzle mass), C_Force (force coefficient). — [same thesis](https://di.univ-blida.dz/jspui/bitstream/123456789/6224/1/32-620-91-1.pdf)

  | M_S | L_d/y* | L/y* | θ* (°) | θ_max (°) | M* | y_S/y* | C_Masse | C_Force |
  |---|---|---|---|---|---|---|---|---|
  | 1.50 | 0.8035 | 2.0143 | 2.4138 | 3.7584 | 1.1515 | 1.0843 | 4.2573 | 0.0615 |
  | 2.00 | 1.3253 | 3.5731 | 5.7466 | 8.2830 | 1.2846 | 1.2984 | 8.6117 | 0.1548 |
  | 3.00 | 2.5457 | 8.3554 | 11.7091 | 15.4652 | 1.4938 | 2.0542 | 29.5532 | 0.2983 |
  | 4.00 | 4.1774 | 16.8144 | 16.1020 | 20.2781 | 1.6427 | 3.2629 | 91.1279 | 0.3818 |
  | 5.00 | 6.3312 | 30.7030 | 19.2722 | 23.5528 | 1.7509 | 4.9748 | 249.8253 | 0.4306 |
  | 6.00 | 9.0873 | 51.9449 | 21.6141 | 25.8751 | 1.8319 | 7.2440 | 611.3943 | 0.4597 |

- **Same tables, HT rows (calorically imperfect air).** Same column order; use for trends only. — [same thesis](https://di.univ-blida.dz/jspui/bitstream/123456789/6224/1/32-620-91-1.pdf)
  - **M_S = 1.5:**
    - T0 = 1000 K: 0.8087, 2.0242, 2.5226, 3.9190, 1.1525, 1.0882, 4.2871, 0.0654
    - T0 = 2000 K: 0.8126, 2.0302, 2.5896, 4.0185, 1.1534, 1.0907, 4.3094, 0.0686
    - T0 = 3000 K: 0.8140, 2.0332, 2.6123, 4.0521, 1.1537, 1.0916, 4.3153, 0.0698
  - **M_S = 2.0:**
    - T0 = 1000 K: 1.3426, 3.6193, 6.0372, 8.6699, 1.2870, 1.3150, 8.8087, 0.1653
    - T0 = 2000 K: 1.3629, 3.6668, 6.2874, 9.0044, 1.2913, 1.3307, 9.0070, 0.1771
    - T0 = 3000 K: 1.3700, 3.6831, 6.3742, 9.1204, 1.2929, 1.3362, 9.0775, 0.1814
  - **M_S = 3.0:**
    - T0 = 1000 K: 2.5958, 8.5759, 12.2434, 16.1048, 1.4949, 2.1145, 31.1362, 0.3162
    - T0 = 2000 K: 2.7300, 9.0508, 13.0807, 17.1019, 1.5114, 2.2351, 34.5225, 0.3490
    - T0 = 3000 K: 2.7817, 9.2351, 13.3874, 17.4638, 1.5174, 2.2818, 35.8782, 0.3625
  - **M_S = 4.0:**
    - T0 = 1000 K: 4.2568, 17.3278, 16.7005, 20.9600, 1.6387, 3.3750, 97.0316, 0.4014
    - T0 = 2000 K: 4.6139, 19.0335, 17.9694, 22.3864, 1.6631, 3.7231, 116.9732, 0.4462
    - T0 = 3000 K: 4.8465, 20.1012, 18.6250, 23.1071, 1.6777, 3.9387, 130.2840, 0.4699
  - **M_S = 5.0:**
    - T0 = 1000 K: 6.4385, 31.6488, 19.8805, 24.2277, 1.7419, 5.1460, 266.3443, 0.4505
    - T0 = 2000 K: 7.0652, 35.4476, 21.3163, 25.7970, 1.7670, 5.7935, 334.9706, 0.4995
    - T0 = 3000 K: 7.6553, 38.7904, 22.2585, 26.7937, 1.7888, 6.3554, 400.9872, 0.5297
  - **M_S = 6.0:**
    - T0 = 1000 K: 9.2344, 53.5639, 22.2279, 26.5444, 1.8190, 7.4928, 652.2294, 0.4799
    - T0 = 2000 K: 10.1522, 60.3605, 23.7059, 28.1349, 1.8419, 8.4865, 831.6022, 0.5302
    - T0 = 3000 K: 11.1976, 67.4405, 24.7883, 29.2554, 1.8670, 9.5064, 1038.9061, 0.5629
- **Zebbiche thesis, Table 5.48: inflection-point position on the axisymmetric MLN wall (x/y*, y/y*).**
  - GP (γ = 1.402):
    - M_S = 1.50: (0.4317, 1.0261)
    - M_S = 2.00: (0.5666, 1.0769)
    - M_S = 3.00: (0.7074, 1.1845)
    - M_S = 4.00: (0.8239, 1.2894)
  - HT, T0 = 1000 K:
    - M_S = 1.50: (0.4333, 1.0274)
    - M_S = 2.00: (0.5703, 1.0811)
    - M_S = 3.00: (0.7137, 1.1945)
  - HT, T0 = 2000 K:
    - M_S = 1.50: (0.4344, 1.0281)
    - M_S = 2.00: (0.5715, 1.0845)
    - M_S = 3.00: (0.7228, 1.2100)
  - HT, T0 = 3000 K:
    - M_S = 1.50: (0.4340, 1.0284)
    - M_S = 2.00: (0.5739, 1.0860)
    - M_S = 3.00: (0.7338, 1.2182)
  - The M_S ≥ 4 HT rows were truncated in my extraction. — [same thesis](https://di.univ-blida.dz/jspui/bitstream/123456789/6224/1/32-620-91-1.pdf)
- **Zebbiche thesis, Table 5.49: M and θ (°) at the inflection point.**
  - GP (γ = 1.402):
    - M_S = 1.50: M = 1.2765, θ = 3.7587
    - M_S = 2.00: M = 1.5032, θ = 8.2892
    - M_S = 3.00: M = 1.8430, θ = 15.4645
    - M_S = 4.00: M = 2.0892, θ = 20.2754
  - HT, T0 = 1000 K:
    - M_S = 1.50: 1.2775, 3.9192
    - M_S = 2.00: 1.5053, 8.6698
    - M_S = 3.00: 1.8423, 16.1035
  - HT, T0 = 2000 K:
    - M_S = 1.50: 1.2785, 4.0187
    - M_S = 2.00: 1.5099, 9.0044
    - M_S = 3.00: 1.8615, 17.0992
  - HT, T0 = 3000 K:
    - M_S = 1.50: 1.2788, 4.0523
    - M_S = 2.00: 1.5120, 9.1199
    - M_S = 3.00: 1.8710, 17.4614
  - The θ at the inflection point equals the θ_max column of Tables 5.42–5.47, to within 0.001°. — [same thesis](https://di.univ-blida.dz/jspui/bitstream/123456789/6224/1/32-620-91-1.pdf)
- **Planar MLN cross-check from the same thesis (Table 5.16, 2D MLN, M_S = 1.50, GP γ = 1.402).** Row values: 1.0365, 2.3512, θ* = 5.9453°, 1.2919, 1.1759, 2.3591, 0.0616. Here θ* = 5.9453° matches ν(1.5)/2. This confirms that the thesis's GP model reproduces the planar MLN result θ_corner = ν_e/2. — [same thesis](https://di.univ-blida.dz/jspui/bitstream/123456789/6224/1/32-620-91-1.pdf)
- **Abada, Kbab & Haif, INCAS Bulletin 16(2) 2024, pp. 3–15** ("Transonic Flow Field Analysis of a Minimum Nozzle Length Rocket Engine", doi:10.13111/2066-8201.2024.16.2.1, CC BY-NC-ND).
  - Table 1, "Numerical values for sizing the nozzles giving the exit Mach number ME", gives axisymmetric MLN results from the HT model. The PDF text extraction scrambled the M* and θ* rows; the de-scrambled values are below.
  - Per M_E:
    - M_E = 2.00: L/y* = 3.6497, M* = 1.2915, θ* = 6.2937°, AE/A* = 1.7711, CMass = 9.0569, CF = 0.1767
    - M_E = 3.00: L/y* = 9.0175, M* = 1.5114, θ* = 13.0812°, AE/A* = 4.9995, CMass = 34.6950, CF = 0.3478
    - M_E = 4.00: L/y* = 19.0236, M* = 1.6632, θ* = 17.9718°, AE/A* = 13.8935, CMass = 117.3712, CF = 0.4440
    - M_E = 5.00: L/y* = 35.4889, M* = 1.7671, θ* = 21.3171°, AE/A* = 33.7028, CMass = 336.1706, CF = 0.4962
  - T0 is not stated next to the table. — [INCAS PDF](https://bulletin.incas.ro/files/abada_kbab_haif__vol_16_iss_2.pdf)
- **Off-the-shelf codes.**
  - MATLAB File Exchange 85253, "Axisymmetric Minimum Length Nozzle Design Tool" (Ludovico Fossà, v1.0.2, Jan 2021), computes the axisymmetric MLN divergent contour for a perfect gas. The page publishes no example numbers and no validation. — [MATLAB FEX 85253](https://www.mathworks.com/matlabcentral/fileexchange/85253-axisymmetric-minimum-length-nozzle-design-tool)
  - A search-result snippet says the tool assumes a normal (straight) sonic line and normalises lengths by throat radius. I did not verify this in the code.

### Inferences
- **The INCAS Table 1 is the HT T0 = 2000 K case.** **[derived]** Its rows match the thesis's HT T0 = 2000 K rows to within about 0.1–2 % (for example, M = 3: 9.0175 vs 9.0508; θ* 13.0812 vs 13.0807; M* 1.5114 identical). It also has AE/A* = 1.7711 at M = 2, where a perfect gas gives 1.6875. Do not use it as a perfect-gas reference.
- **GP exit radii check out against 1-D A/A*.** **[derived]** (y_S/y*)² = 1.6858 at M = 2 (isentropic γ = 1.4: 1.6875) and 4.2197 at M = 3 (γ = 1.4: 4.2346; γ = 1.402 is slightly lower). So y* is the throat radius and the sonic line is straight (uniform M = 1 at the throat). This makes rockmoc's existing A/A* check equivalent to this column.
- **Suggested rockmoc test.** With γ = 1.402 and a straight sonic line, compare L/y*, θ*, θ_max and y_S/y* against the GP row of the Zebbiche table. Use a tolerance of about 1–2 % on L/y*: the thesis's Tables 5.35–5.41 show the results depend on mesh and Kernel refinement, and the INCAS rerun differs from the thesis by about 0.5 % for the same model.
- **Corner angle vs maximum wall angle.** In this axisymmetric MLN formulation the corner angle θ* is smaller than the maximum wall angle θ_max. The wall keeps turning outward after the corner and reaches θ_max at an inflection point near x/y* ≈ 0.4–0.8. rockmoc should report and test both values. A test that assumes "θ_max at the sharp corner" will disagree with this source. The task brief's phrase "θ_max at the sharp corner" corresponds to Zebbiche's θ*. Whether Argrow & Emanuel's tabulated wall angle is the corner angle or the maximum is unverified.

### Gaps
- No Argrow & Emanuel (1988) numbers were found in open sources. The ASME paper is paywalled. The 2025 Advances in Aerodynamics paper cites it without reproducing it. Argrow's Oklahoma thesis was not located online.
- Zebbiche & Youbi / Zebbiche JBIS 63 (2010) 171–191 ("Supersonic axisymmetric MLN conception at high temperature with application for air") was not accessible; only the bis-space.com product page exists. Its tables probably match the thesis.
- **Unsourced length table in evand/nozzle.** The repo's references README has a "Nozzle Length Reference Data" table with "L_MLN/r* (approx)" values of about 1.5 / 3.5 / 8 / 18 for M = 2 / 3 / 4 / 5. It gives no source, and the values are about half of Zebbiche's L/y* (3.57 / 8.36 / 16.8 / 30.7). Treat it as unreliable, possibly a planar half-height or a different normalisation. — [evand/nozzle references README](https://github.com/evand/nozzle/blob/main/docs/references/README.md)
- I did not run or inspect the milo-thatch repo or the MATLAB FEX 85253 code, so neither has verified printed numbers.

## Q2. Theoretical facts: θ_max vs ν_e, and axisymmetric vs planar length

### Takeaway
For a planar MLN with a straight sonic line, θ_w,max = ν_e/2 exactly. For an axisymmetric MLN there is no algebraic relation, and the corner angle θ* is well below ν_e/2. In Zebbiche's GP tables θ* is roughly 0.18–0.25 × ν_e and θ_max roughly 0.25–0.32 × ν_e **[derived ratios]**. Axisymmetric MLNs are also much shorter, relative to throat radius, than planar MLNs relative to throat half-height.

### Cited Findings
- Planar MLN: θ_max = ν(M_e)/2. This is used as a test in the evand/nozzle code (Anderson, Modern Compressible Flow, Eq. 11.33; γ = 1.4, M_e = 2.4 gives θ_max ≈ 18.38°). — [evand/nozzle VALIDATION.md](https://github.com/evand/nozzle/blob/main/docs/VALIDATION.md)
- A secondary source (search-result snippet of a ResearchGate paper; full text returned HTTP 403) states that for axisymmetric nozzles the final wall angle has no algebraic relationship, because the characteristic lines vary through the flow field. It also states that a straight sonic line assumption gives shorter nozzles than other MOC approaches. — [ResearchGate: Short nozzles design for real gas supersonic flow using MOC](https://www.researchgate.net/publication/358224507_Short_nozzles_design_for_real_gas_supersonic_flow_using_the_method_of_characteristics) (unverified: snippet only)
- The Zebbiche thesis states that the axisymmetric MLN wall has an inflection point and the 2D MLN wall does not. — [Zebbiche thesis](https://di.univ-blida.dz/jspui/bitstream/123456789/6224/1/32-620-91-1.pdf)
- evand/nozzle's design notes say the simplification θ_w = θ_I is "Only exact for 2D planar; approximate for axisymmetric". They also say the axisymmetric source term Q = sinθ·sinμ/(y·cos(θ±μ)) makes the wall angle diverge in a naive wall-point algorithm, and that there is "No verified reference for the exact wall point algorithm in axisymmetric flow". — [evand/nozzle axisymmetric_mln.md](https://github.com/evand/nozzle/blob/main/docs/algorithms/axisymmetric_mln.md)

### Inferences
- **Corner and maximum angles vs ν_e/2.** **[derived]** The ν_e values below use γ = 1.4; γ = 1.402 changes them by less than 0.1°.

  | M_e | ν_e (°) | ν_e/2 (°) | Zebbiche GP θ* (°) | Zebbiche GP θ_max (°) |
  |---|---|---|---|---|
  | 2 | 26.38 | 13.19 | 5.75 | 8.28 |
  | 3 | 49.76 | 24.88 | 11.71 | 15.47 |
  | 4 | 65.79 | 32.89 | 16.10 | 20.28 |

  Both axisymmetric angles are well under ν_e/2.
- **Corner Mach number.** M* is consistent with ν(M*) = θ*. For example, ν(1.4938) ≈ 11.7° at M_S = 3. So the corner fan in this formulation turns the flow from 0 to θ*, exactly as in planar flow.
- **Length.** For comparison, the planar MLN length in throat half-heights at M = 2 from the thesis's 2D Table 5.17 was not transcribed (only the M = 1.5 2D row was). The axisymmetric L/y* = 2.0143 is shorter than the planar 2.3512 at M_S = 1.5, consistent with the expected shorter axisymmetric MLN.

### Gaps
- I found no closed-form or published asymptotic relation for θ*(M_e, γ) in axisymmetric MLNs.
- Values for γ other than 1.402 for perfect-gas axisymmetric MLNs were not found in open sources. The Argrow & Emanuel paper reportedly covers other γ (unverified).

## Q3. Published wall-coordinate tables for axisymmetric MOC (ideal / wind-tunnel) nozzles

### Takeaway
The canonical source is Sivells, AEDC-TR-78-63 (Dec 1978), DTIC AD-A062944. Its sample is a Mach 4 axisymmetric nozzle: γ = 1.4, η = 8.67°, RC = 6, inviscid exit radius 12.25 in. The scanned report's OCR text is poor. The MIT-licensed typeset re-creation (aldorona/contur) includes the report PDF and a clean `docs/output.txt` with full wall-coordinate tables. However, the evidence indicates that this output file is a planar run (JD = −1, BMACH = 3.2), not the report's axisymmetric sample, and evand/nozzle tests quote it as the "Sivells Mach 4" case. The 2025 Advances in Aerodynamics paper gives dimensional axisymmetric hypersonic nozzle sizes (r* = 0.05 m) for M = 6, 8, 10 and γ = 1.4, 1.34.

### Cited Findings
- **Sivells source and DTIC location.** Sivells JC, "A computer program for the aerodynamic design of axisymmetric and planar nozzles for supersonic and hypersonic wind tunnels", AEDC-TR-78-63, ARO Inc., Dec 1978. The scanned listing is at https://apps.dtic.mil/dtic/tr/fulltext/u2/a062944.pdf (DTIC returned "request blocked"/403 to my fetches). The aldorona/contur repo contains the F77 listing (OCR-typeset by F.L. Zavalan and A. Rona), an F90 port, `sivells/Sivells.pdf` (3.8 MB), and `docs/output.txt`. The README claims that "run in May 2019 on this test case and the same output was obtained". — [aldorona/contur](https://github.com/aldorona/contur)
- **Sivells report §7.0, sample design.** "The design of a Mach 4 axisymmetric nozzle is selected to illustrate use of the computer program… An axisymmetric nozzle is specified by leaving JD blank (JD = 0)… The computed value of FMACH is 3.0821543, which is greater than the value of BMACH specified on Card 3; therefore, BMACH also becomes 3.0821543. The negative value of SF means that the inviscid exit radius of the nozzle is 12.25 in. The value of PP means that the inflection point will be 60 in. downstream of an arbitrary point."
  - Input cards: γ = 1.4, R = 1716.563, ETAD = 8.67, RC = 6., BMACH = 3., CMC = 4., SF = −12.25, PP = 60.
  - Readable sample-output values from the report's downstream-contour page (OCR, partly garbled): SCALE FACTOR (SF) = 21.750386; BMACH = 3.08215; AMACH = 2.2878437; XA = 1.4559664; XB = 2.1398501; XC = 6.3776623; XB(IN) = 76.9263852; XC(IN) = 132.3126057; XD(IN) = 179.7566517; YD(IN) = 12.2500000. — [Sivells.pdf in aldorona/contur](https://github.com/aldorona/contur/blob/master/sivells/Sivells.pdf)
- **Sivells report, transonic note.** "the second derivative of velocity ratio at the sonic point is negative for values of R less than 11.767 for planar flow and 10.525 for axisymmetric flow." Here R is the throat radius of curvature ratio. — [Sivells.pdf](https://github.com/aldorona/contur/blob/master/sivells/Sivells.pdf)
- **aldorona/contur `docs/output.txt` header values** (all copied from the file):
  - "GAMMA= 1.4000 INFLECTION ANG. (ETA)= 8.6700 DEGREES RAD. OF CURV. (RC)= 6.000000 SCALE FACTOR (SF)= 7.55258257"
  - "Y*=0.15132005 RMASS=0.99937483 … EMACH= 1.44128 FMACH= 2.0417400"
  - "AMACH= 2.7550181 XA= 3.3152830 XB= 5.1209575 XBC= 13.1696977 XC= 18.2906552 XD= 24.5724861"
  - "XA(IN)= 60.0000000, YA(IN)= 3.8180812, XB(IN)= 73.6375059, XC(IN)= 173.1027355, XD(IN)= 220.5467815, YD(IN)= 12.2500000"
  - Mach 1 at 41.8494756 in.

  Selected upstream inviscid wall points (X, Y nondimensional; M; flow angle °):
  - pt 1: 0.89662098, 0.15141471, 1.0639748, 0.0
  - pt 11: 0.95415029, 0.15319460, 1.1959967, 3.4577595
  - pt 21: 1.0321348, 0.16024778, 1.3361384, 6.4690204
  - pt 31: 1.1666232, 0.17818468, 1.5185491, 8.3092745
  - pt 41: 1.3563968, 0.20683110, 1.7356275, 8.6700000

  Downstream:
  - pt 41: 3.3152830, 0.50553320, M 2.7550181, 8.67°
  - pt 149 (exit): 24.572486, 1.6219617, M 4.0000000, 0.0°

  Boundary-layer-corrected contour: length = 178.8245270 in; Y(IN) = 15.113692 at X = 172 in. — [docs/output.txt](https://github.com/aldorona/contur/blob/master/docs/output.txt)
- **evand/nozzle uses this output as its "Sivells Mach 4" axisymmetric validation.** Its test file says "y* = 0.15132005 (y*/y0, throat half-height / reference height)". It also notes that CONTUR's η = 8.67° is a design parameter, not ν/2. — [evand/nozzle tests/test_contur_equivalence.py](https://github.com/evand/nozzle/blob/main/tests/test_contur_equivalence.py)
- **Advances in Aerodynamics 2025** (s42774-025-00244-1), "A modular method of characteristics framework for controllable-length axisymmetric high-speed wind tunnel nozzle design". Code: https://github.com/PERRY0125/Hypersonic-nozzle-framework (v1.0.0).
  - Table 1, "Nozzle dimensions for different Mach numbers and γ values", for throat radius 0.05 m, "maximum expansion angles of 11° and 12°", P0 = 15 MPa, T0 = 1200 K:

    | Mach | γ | Initial expansion length (m) | Total length (m) | Exit radius (m) |
    |---|---|---|---|---|
    | 8 | 1.4 | 3.3559 | 8.8362 | 0.6901 |
    | 8 | 1.34 | 4.7469 | 12.0274 | 0.9172 |
    | 10 | 1.4 | 5.2140 | 16.7145 | 1.1563 |
    | 10 | 1.34 | 8.2645 | 25.3466 | 1.7204 |

  - Mach 6 baseline: throat radius 0.05 m, throat arc radius 1.045 m, maximum expansion angle 10°, total length 3.999 m, outlet radius 0.364 m. One later passage says 0.464 m.
  - The downstream "damping section" references the MLN design method, and its wall tracing uses "Zebbiche's method". — [SpringerLink article](https://link.springer.com/article/10.1186/s42774-025-00244-1); [Table 1](https://link.springer.com/article/10.1186/s42774-025-00244-1/tables/1)
- **Other primary sources with tabulated coordinates** (titles only; their tables were not retrieved), as cited in the 2025 paper's references:
  - Cresci RJ (1958), "Tabulation of coordinates for hypersonic axisymmetric nozzle: Part 1, analysis and coordinates for test section Mach numbers of 8, 12 and 20", WADC TN 58-30
  - Crown JC (1948), NACA TN-1651
  - Sivells (1970), J. Spacecraft Rockets 7(11), with hypersonic validation at M = 8, 10, 20 (per evand's README)
  - Guentert & Neumann, NASA TR-R-33 (1959)
  - Beckwith & Moore, NACA TN-3322

  — [Advances in Aerodynamics reference list](https://link.springer.com/article/10.1186/s42774-025-00244-1); [evand/nozzle references](https://github.com/evand/nozzle/blob/main/docs/references/README.md)

### Inferences
- **`docs/output.txt` appears to be a planar run.** **[derived]** It differs from the report's axisymmetric sample in SF (7.55258257 vs 21.750386), FMACH/BMACH (2.0417/3.2 vs 3.0821543/3.0821543) and XD(IN) (220.55 vs 179.76 in). Its exit-to-throat height ratio, 1.6219617/0.15141471 = 10.712, equals the 1-D A/A*(M = 4, γ = 1.4) = 10.719, which is the planar relation. For an axisymmetric nozzle the radius ratio would be √10.72 ≈ 3.27. The contur README says that setting "JD=-1 and BMACH=3.2" gives a planar contour, which matches the values in this output. **So `docs/output.txt` (and evand's "Sivells Mach 4" test values derived from it) appear to be the PLANAR variant.** It remains a good planar ideal-contour reference. Do not use it as an axisymmetric one.
- **Axisymmetric Sivells sample: derived throat size.** **[derived, unverified]** The inviscid exit radius is 12.25 in, and SF = 21.750386 gives nondimensional y_exit = 0.5632. The throat radius should be about 12.25/3.27 ≈ 3.74 in, slightly adjusted by the RMASS mass-deficit factor. Confirm this against a clean copy of the report's sample output before using it.
- **Advances in Aerodynamics exit radii: derived checks.** **[derived]** All values are consistent with 1-D isentropic A/A*:
  - M8, γ1.4: (0.6901/0.05)² = 190.5 vs A/A* = 190.1
  - M10, γ1.4: (1.1563/0.05)² = 534.8 vs 535.9
  - M6: 0.364/0.05 = 7.28, (7.28)² = 53.0 vs 53.18
  - So 0.364 m is correct and "0.464 m" is a typo.
  - Nondimensional lengths: L/r* = 176.7 (M8, γ1.4), 334.3 (M10, γ1.4), 80.0 (M6).
  - Which expansion angle (11° or 12°) belongs to which Mach is not stated explicitly.
  - These are ideal (not minimum-length) contours, useful for the ideal-contour designer only at the A/A* level unless the PERRY0125 code is run.

### Gaps
- No clean, machine-readable wall-coordinate table for Sivells' axisymmetric Mach 4 sample was obtained. The OCR of `Sivells.pdf` (pp. ~70+) is garbled, and DTIC blocked direct download.
- Next steps to close this:
  - Rerun the aldorona F90 code with the report's exact Card 3 (BMACH = 3., JD = 0).
  - Confirm that SF = 21.750386 and XD(IN) = 179.7566517 are reproduced, then use that output.
- Cresci WADC TN 58-30, NASA TR-R-33 and Beckwith & Moore NACA TN-3322 coordinate tables were not fetched in this session. NTRS hosts the NACA/NASA items; locations unverified.
