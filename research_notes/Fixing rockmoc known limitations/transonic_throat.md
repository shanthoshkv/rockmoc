# Transonic throat-flow solutions for tight axisymmetric throats (Rc/Rt < 1.5)

Notation: R = Rc/Rt, lengths in units of Rt, x measured from the throat plane, velocities u, v = V/a*. Several equations below were rebuilt from OCR'd scans and checked numerically in Python; each item says whether it was verified from a source or reconstructed. Scratch scripts: `C:\Users\kvsha\AppData\Local\Temp\claude\C--Users-kvsha-Documents-Learning-rockets-rockmoc\c646cc12-ba86-4116-b614-06e8f6f6486a\scratchpad\kl.py` (series) and `...\scratchpad\res.py` (residual check).

## Q1. Kliegel & Levine (1969): the series, coordinates, sonic line and Cd

### Takeaway
Kliegel-Levine (KL) did not derive new coefficients. They took Hall's 3rd-order cylindrical-coordinate series in 1/R (with corrections to Hall) and **re-expanded it in ε = 1/(R+1)**. Hall's u1, u2, u3, v1, v2, v3 are reused unchanged, combined with binomial weights. Levine & Coats (1970) later showed the KL series does not satisfy the equations exactly, so it is "an empiricism". Even so, it matches tight-throat data down to R ≈ 0.4. It is a drop-in upgrade for rockmoc's `hall_uv`: about 30 lines, no new solver.

### Cited Findings
- **Where the toroidal idea came from.** Kliegel & Levine (Ref. 10 in the report below) argued that the R ≳ 1 limit came from the cylindrical coordinates, and proposed both toroidal coordinates and the parameter 1/(1+R), which never exceeds 1. "Hall's solution was recast as series in 1/(1+R) in Reference 10." The results were "very encouraging" against Cuffel et al. data at R = 0.625. — [Levine & Coats, "Transonic Flow in a Converging-Diverging Nozzle", Dynamic Science SN-168, NASA contract NAS7-756, Sept 1970 (NTRS 19710000970), PDF](https://ntrs.nasa.gov/api/citations/19710000970/downloads/19710000970.pdf)
- **Levine & Coats verdict on KL.** "The results obtained using the method of Reference 10 show excellent agreement with the data for R down to approximately .4; however, this must be considered fortuitous." The KL series is not the toroidal solution transformed back to cylindrical coordinates, and it "do[es] not satisfy the differential equations of motion in cylindrical coordinates … actually an empiricism which agrees quite well with the data … it still represents a useful and unique engineering tool". Extending KL to 4th order made it worse, and 5th order diverged. — [same report, pp. 39–41](https://ntrs.nasa.gov/api/citations/19710000970/downloads/19710000970.pdf)
- **True toroidal-coordinate solutions (Levine & Coats) did not help.** They solved recursively to 5th order with ε = (1+2R)^½/(1+R) or ε = (1/(1+R))^½, including variable γ. "None of them appear to be capable of yielding valid, accurate solutions for small R." Hall's cylindrical 3rd order diverges below R ≈ 1.5, and all higher orders diverge below R ≈ 2. The explanation offered: for small R the upstream geometry affects the throat flow, so any expansion that uses only the local throat geometry fails. — [same report, Abstract and pp. 35–41](https://ntrs.nasa.gov/api/citations/19710000970/downloads/19710000970.pdf)
- **Toroidal transformation** (from the same report, equations 1–3; OCR-damaged, so I re-derived it — see Inferences):
  - r/Rt = (1+2R)^½ sinh η / (cos ξ + cosh η)
  - z/Rt = (1+2R)^½ sin ξ / (cos ξ + cosh η)
  - The throat wall is the line η = η_w.
  - Velocity conversion: v_r = [u sinh η sin ξ + v (cosh η cos ξ + 1)]/(cos ξ + cosh η), with v_z analogous. **Not verified**: the scan is illegible here.
- **KL series as transcribed from the Blida master's thesis** (eqs. II.31–II.38), which quotes Kliegel & Levine. Uses Hall's z = x·√(2R/(γ+1)), r = y/yt, ε = 1/(R+1):
  - u = 1 + ε u1 + ε² (u1 + u2) + ε³ (u1 + 2u2 + u3)
  - v = √((γ+1)ε/2) · [ ε v1 + ε² (3/2 v1 + v2) + ε³ (15/8 v1 + 5/2 v2 + v3) ]
  - — [Kbab K., "Développement de méthodes de conception de tuyères supersoniques de propulsion", Master's thesis, Univ. Saad Dahleb Blida, 2012, pp. 31–33, PDF](https://di.univ-blida.dz/jspui/bitstream/123456789/3283/1/036-2012.pdf)
  - These weights are exactly the binomial re-expansion of Hall: 1/R = ε/(1−ε) and (1/R)^{3/2} = ε^{3/2}(1 + 3/2 ε + 15/8 ε² …). That confirms the KL construction.
- **Hall / KL coefficients**, same thesis. Signs were lost in the OCR and restored from rockmoc's existing Hall implementation and from the irrotationality check:
  - u1 = r²/2 − 1/4 + z
  - v1 = r³/4 − r/4 + r z
  - u2 = (2γ+9)/24 r⁴ − (4γ+15)/24 r² + (10γ+57)/288 + z(r² − 5/8) − (2γ−3)/6 z²
  - v2 = (γ+3)/9 r⁵ − (20γ+63)/96 r³ + (28γ+93)/288 r + z[(2γ+9)/6 r³ − (4γ+15)/12 r] + r z²
  - u3 = (556γ²+1737γ+3069)/10368 r⁶ − (388γ²+1161γ+1881)/2304 r⁴ + (304γ²+831γ+1242)/1728 r² − (2708γ²+7839γ+14211)/82944
    + z[(52γ²+51γ+327)/384 r⁴ − (52γ²+75γ+279)/192 r² + (92γ²+180γ+639)/1152]
    + z²[−(7γ−3)/8 r² + (13γ−27)/144] + (4γ²−57γ+27)/144 z³
  - v3 = (6836γ²+23031γ+30627)/82944 r⁷ − (3380γ²+11391γ+15291)/13824 r⁵ + (3424γ²+11271γ+15228)/13824 r³ − (7100γ²+22311γ+30249)/82944 r
    + z[(556γ²+1737γ+3069)/1728 r⁵ − (388γ²+1161γ+1881)/576 r³ + (304γ²+831γ+1242)/864 r]
    + z²[(52γ²+51γ+327)/192 r³ − (52γ²+75γ+279)/192 r] − (7γ−3)/12 r z³
  - — [Blida thesis, eqs. II.36–II.38](https://di.univ-blida.dz/jspui/bitstream/123456789/3283/1/036-2012.pdf)
  - u1, u2, v1, v2 match `src/rockmoc/moc/transonic.py::hall_uv` exactly.
  - The OCR printed "1735" in the z·r⁵ term of v3. I used 1737, the value irrotationality requires.
- **KL also corrected an error in Hall's Cd formula** (about 0.06 % effect) and gave a 1/(1+R) Cd series that "converges for all values of R". — [Li et al., FLOMEKO 2019 (IMEKO TC9-2019-028), PDF](https://imeko.org/publications/tc9-2019/IMEKO-TC9-2019-028.pdf). The formulas are under Q5.
- **Sonic line.** A secondary source says the KL sonic line is parabolic, starts at the throat, and that the 1/(R+1) expansion moves it downstream relative to Hall. — [Abada/Kbab et al., "Application of De-Laval nozzle transonic flow field computation approaches" (academia.edu; full text returned 403, only the snippet was read)](https://www.academia.edu/6516384/Application_of_De_Laval_Nozzle_Transonic_Flow_Field_Computation_Approaches)

### Inferences
- **My Python check, γ = 1.4, throat plane (x = 0):**

  | R | Method | Axis M | Wall M |
  |---|---|---|---|
  | 0.625 | Hall 2nd order | 1.300 | 1.114 |
  | 0.625 | Hall 3rd order | 0.252 | 2.20 |
  | 0.625 | **KL 3rd order** | **0.790** | **1.314** |
  | 1.0 | KL | 0.837 | 1.226 |
  | 1.0 | Hall 2nd order | 0.996 | — |
  | 2.0 | KL | 0.897 | 1.129 |
  | 2.0 | Hall 2nd order | 0.925 | — |
  | 2.0 | Hall 3rd order | 0.873 | — |

  - At R = 0.625, KL's axis Mach is close to the measured ≈ 0.8 (Cuffel). Its wall Mach is about 7 % below the quoted 1.4.
  - The fact that the axis value reproduces the measurement is the strongest evidence that my u3 constant term and the KL form are transcribed correctly.
- **Consistency checks:**
  - With the corrected v3 z³ term, the 3rd-order series satisfies irrotationality (∂u/∂r = ∂v/∂x) to machine precision.
  - The continuity residual drops from ε^1.3 (Hall 2nd order) to ε^2.6 (Hall 3rd order) on my test patch. This test patch differs from rockmoc's `_potential_residual`.
- **Toroidal parameters, re-derived by me rather than read from the scan.**
  - The focal radius is a = Rt√(1+2R).
  - At the wall, cosh η_w = (1+R)/R and sinh η_w = √(1+2R)/R, so η_w = ln[(1+R+√(1+2R))/R].
  - Check: at ξ = 0, r = a·tanh(η_w/2) = Rt.
  - Only needed if a true toroidal solver is ever written. Not needed for KL.
- **Axis sonic point from KL:** x*/Rt ≈ +0.173 (R = 0.625), +0.161 (R = 1), +0.149 (R = 1.5), +0.139 (R = 2).
- **KL wall tangency** on the arc is within 1° for −0.1 ≤ x ≤ 0.1 at R = 0.625, and within 0.15° at R = 1.
- **Caution for the start line.** At R = 0.625 the KL wall Mach rises steeply downstream: 1.45 at x = 0.05, 1.62 at 0.10, 2.17 at 0.20. A vertical start line at the axial Mach 1.02 station (x ≈ 0.2) would therefore carry a wall Mach that is probably too high.
  - Best option: start on a non-vertical line, for example from the axis point at M ≈ 1.01 back to the wall near x ≈ 0–0.05. KL is most trustworthy close to the throat plane.
  - Simpler option: start from the throat-plane-to-sonic region and let MOC march.
  - Either way, validate against Cuffel's wall data. This is my inference, not tested against data.
- **Implementation is tiny.**
  - `kl_uv` = `hall_uv` with ε = 1/(R+1), plus u3 and v3, plus the binomial weights. Keep z = x√(2R/(γ+1)): in the KL form only the outer parameter changes.
  - Unverified risk: the thesis OCR scaling line could also be read as √(2(R+1)/(γ+1)). The axis and wall values on the throat plane do not depend on this choice, but the axial gradients do. Check it against the original paper or against Cuffel's axial distribution.

### Gaps
- I could not access the original KL paper (AIAA J 7(7):1375–1378, [doi 10.2514/3.5355](https://arc.aiaa.org/doi/10.2514/3.5355)). That leaves unverified:
  - the exact z-scaling;
  - whether KL kept Hall's u3/v3 including Hall's original error, or a corrected form;
  - the sonic-line formula KL gives.
- Zucrow & Hoffman Vol. 2 is not online. I could not confirm whether it reproduces KL.

## Q2. Dutton & Addy

### Takeaway
I found only metadata. Dutton & Addy's series is still an inverse-powers-of-R throat expansion (axisymmetric, planar and annular nozzles), implemented in the FORTRAN code TRANNOZ. Nothing found shows it fixes R < 1. The 1970 analysis above argues that no local expansion can.

### Cited Findings
- **TRANNOZ.** "A computer program for analysis of transonic throat flow"; Univ. Illinois report UILU-ENG 80-4005, April 1980; evaluates the Dutton–Addy series for axisymmetric, planar and annular nozzles. — [Illinois Experts entry](https://experts.illinois.edu/en/publications/trannoz-a-computer-program-for-analysis-of-transonic-throat-flow-/)
- **Annular-nozzle paper.** Dutton & Addy, "Transonic flow in the throat region of annular supersonic nozzles", AIAA J 20(9):1236–1243, 1982. Agreement with pressure measurements is "quite good through a significant region of the throat", with discrepancies near the walls in the throat inlet region. — [Illinois Experts entry](https://experts.illinois.edu/en/publications/transonic-flow-in-the-throat-region-of-annular-supersonic-nozzles/)

### Inferences
- The axisymmetric Dutton–Addy 1981 paper (AIAA J 19(6)) is very likely another higher-order 1/R expansion, so it is not a remedy for R = 0.625. **Not verified.**

### Gaps
- Equations, validity range and the Cuffel comparison for Dutton & Addy (1981) were not found. The paper is AIAA-paywalled, and the TRANNOZ report was not found online.

## Q3. Other approaches and which is simplest

### Takeaway
For R < 1, the only rigorously valid route is a numerical transonic solver: a time-marching Euler solution, or the methods of Prozan, Serra, or Moretti. The simplest acceptable route for rockmoc is the **KL 3rd-order series**:
- use it for 0.4 ≲ R < 1.5 and keep Hall 2nd order for R ≥ 1.5;
- below R ≈ 0.4, add an optional small 2-D axisymmetric time-marching (MacCormack) solver over the throat region and extract the start line from it.

### Cited Findings
- **Prozan's finite-difference solver.** Prozan (Lockheed Huntsville, NAS7-743) used finite-difference (time-dependent) techniques on eight JPL contours with steep inlets (75°, 90°) and **R = 0.5 and 0.25**. "The finite-difference techniques used in the analysis performed very well, particularly when the severity of the contours is considered." The JPL technical manager was R. F. Cuffel. — [Prozan, "Transonic Flow in a Converging-Diverging Nozzle", Final Report (NTRS 19700017856), PDF](https://ntrs.nasa.gov/api/citations/19700017856/downloads/19700017856.pdf)
- **Numerical methods in 1970.** They "are not subject to the limitation on R", but at that time were expensive and subject to instability and to iteration on the choked mass-flow singularity. — [Levine & Coats 1970, Introduction](https://ntrs.nasa.gov/api/citations/19710000970/downloads/19710000970.pdf)
- **Levine & Coats Section IV.** They proposed generating a *subsonic* start line from the local transonic expansion and marching the subsonic field, to avoid iterating on the mass-flow singularity. They restrict this to R > 1. — [same report, Section IV/V](https://ntrs.nasa.gov/api/citations/19710000970/downloads/19710000970.pdf)
- **Kliegel & Quan.** "Convergent-Divergent Nozzle Flows" is an expansion approach cited as Ref. 6 in the same report. — [Levine & Coats reference list](https://ntrs.nasa.gov/api/citations/19710000970/downloads/19710000970.pdf)
- **TDK.** NASA's TDK code is described as a 2-D MOC with finite-rate kinetics and a "Transonic module". I did not find documentation of which transonic start it uses. — [NTRS 19960029260](https://ntrs.nasa.gov/citations/19960029260); TDK engineering manual NTRS 19860007470 (the PDF download failed)
- **Series methods generally.** They are effective for Rc/Rt > 1 and unsuitable below 1. — [Brazilian JATM study citing Wear (1971), via search snippet](https://scielo.br/j/jatm/a/JnYTP8jWqdVVWyQJdMjhRMS/?lang=en)
- **Sonic-line shape for small Rc.** For small radius of curvature, the sonic line starts upstream of the throat at the wall and ends downstream on the axis. — [IIETA paper, via search snippet](https://iieta.org/download/file/fid/4245)

### Inferences
- **Ranking for rockmoc:**
  1. KL series: hours of work, closed form.
  2. Pseudo-1-D axis plus series: no gain.
  3. A 2-D axisymmetric MacCormack time-marching Euler solver on a small structured grid from about x = −1.5 Rt to +0.5 Rt: about 300 lines of numpy. Its upstream boundary is p0/T0 at the inlet; at the outlet it is supersonic extrapolation. It is robust for any R and any inlet angle, at the cost of seconds of runtime.
  4. A full streamline-curvature method: more complex, not recommended.
- Since KL is "fortuitous" and its error grows as R → 0, gate it with a warning below R ≈ 0.4–0.5. That is the limit where Levine & Coats say agreement holds.

### Gaps
- No open-source Python implementation of KL or Dutton–Addy was found. The GitHub MCP was unavailable this session, and web search surfaced none.
- No description was found of how TDK or OpenFOAM practice sets start lines.

## Q4. Validation data

### Takeaway
The primary tight-throat dataset is Cuffel, Back & Massier (1969), R = 0.625. Back, Massier & Gier (1965) covers R = 2. Levine & Coats Fig. 4 plots throat wall velocity against 1/(1+R) from four datasets. I found no machine-readable tables online.

### Cited Findings
- **Datasets in Levine & Coats Fig. 4** (3rd-order throat wall velocity against 1/(1+R)): Shelton; Cuffel (private communication); Cuffel, Back & Massier (AIAA J 1969); Smith & Matz. — [Levine & Coats 1970, Fig. 4](https://ntrs.nasa.gov/api/citations/19710000970/downloads/19710000970.pdf)
  - The OCR'd axis tick labels run 1.0–2.0 in u_w and 0–1.0 in 1/(1+R). The data points themselves could not be digitised from the OCR text.
- **Other papers citing Cuffel et al.** for transonic Mach distribution and sonic line: [Tolentino et al., FME Transactions 50(2), 2022](https://www.mas.bg.ac.rs/_media/istrazivanje/fme/vol50/2/5_s.l.tolentino_et_al.pdf).
- **rockmoc's existing R = 2 throat wall-pressure taps** (from NASA TN D-5377) are p/p0 = 0.5733, 0.46, 0.3733, 0.2785 at x/Rt = −0.117, 0, 0.087, 0.171. — `website/make_site_data.py:114` in this repo

### Inferences
- **Values for an R = 0.625 regression test**, taken from the task brief (not re-verified here): axis M ≈ 0.8 and wall M ≈ 1.4 at the throat plane, Cd = 0.985. KL gives 0.790 / 1.314 / 0.982 (see Q5). A test should be tolerance-based.

### Gaps
- I could not retrieve the tabulated Cuffel 1969 or Back, Massier & Gier 1965 data:
  - the AIAA J papers are paywalled;
  - the related JPL reports were not located on NTRS this session;
  - the Levine & Coats PDF is a scanned bitmap with no text layer (the text used here came from OCR), so the figure data could not be digitised.

## Q5. Discharge-coefficient correlations (inviscid, function of R and γ)

### Takeaway
There are three closed-form inviscid Cd series, as transcribed in FLOMEKO 2019 (R = Rc/Rt = 2Rc/d). KL's version stays accurate for small R. At R = 0.625 it gives Cd = 0.9817, against a measured 0.985.

### Cited Findings
Source for all three formulas: [Li et al., FLOMEKO 2019, eqs. 1–3](https://imeko.org/publications/tc9-2019/IMEKO-TC9-2019-028.pdf). Signs were partly reconstructed from the OCR layout and checked numerically (see Inferences).

- **Hall (1962), original:**
  1 − Cd = (γ+1)/R² · [ 1/96 − (8γ+21)/(4608 R) + (754γ²+1971γ+2007)/(552960 R²) ]
- **Hall, as corrected by Kliegel–Levine:**
  1 − Cd = (γ+1)/R² · [ 1/96 − (8γ+21)/(2304 R) + (754γ²+2123γ+2553)/(552960 R²) ]
- **Kliegel–Levine, series in 1/(1+R)**, which "converges for all values of R":
  1 − Cd = (γ+1)/(1+R)² · [ 1/96 − (8γ−27)/(2304 (1+R)) + (754γ²−757γ+3633)/(276480 (1+R)²) ]
- **Hall Cd for small R.** Hall's Cd becomes non-physical (negative) for sufficiently small R. — [same source](https://imeko.org/publications/tc9-2019/IMEKO-TC9-2019-028.pdf)
- **Cross-check of Hall original.** It gives Cd2 = 0.99859 for γ = 1.4 at R = 4 (Rc = 2d); my evaluation of the original Hall formula reproduces 0.99859. — [Li et al., FLOMEKO 2013 (IMEKO TC9-2013-105)](https://imeko.org/publications/tc9-2013/IMEKO-TC9-2013-105.pdf)
- **Viscous correction.** Total Cd is commonly split as Cd = Cd1 (viscous boundary layer, ~b/√Re; Geropp laminar solution) × Cd2 (inviscid, from the formulas above). PTB uses Kliegel's solution for the geometric parameter a. — [IMEKO TC9-2013-105](https://imeko.org/publications/tc9-2013/IMEKO-TC9-2013-105.pdf); [Cao/Mickan PTB, IMEKO TC9-2019-038](https://www.imeko.info/publications/tc9-2019/IMEKO-TC9-2019-038.pdf) (redirected and not read)
- **Tang & Fenn.** "Experimental Determination of the Discharge Coefficients for Critical Flow through an Axisymmetric Nozzle", AIAA J 16(1):41–46, 1978. Only the citation was found, no correlation. — [search snippet via Springer paper s12206-010-0627-6](https://link.springer.com/article/10.1007/s12206-010-0627-6)

### Inferences
- **My evaluation, γ = 1.4:**

  | R | KL Cd formula | Measured | KL velocity field integrated over throat plane | Hall-corrected formula |
  |---|---|---|---|---|
  | 0.625 | 0.9817 | 0.985 (Cuffel) | 0.9812 | — |
  | 1.0 | — | — | 0.9894 | — |
  | 2.0 | — | — | 0.9962 | 0.9960 |

  The agreement between the KL formula and my independently integrated KL velocity field (to 5e-4 at R = 0.625 and 2e-4 at R = 2) cross-validates both the Cd signs and the u/v series transcription.
- **For rockmoc:** use the KL Cd formula as the inviscid Cd for any R, and optionally multiply by a Geropp-type viscous term.

### Gaps
- Gas-dynamic Cd tables against R from Tang & Fenn and from Cuffel were not retrieved.
- The Geropp formula's exact coefficients were not verified (the OCR was garbled).
