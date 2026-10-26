# Closing rockmoc's six physics gaps, ranked

Four of rockmoc's six known limitations can be fixed with closed-form or low-cost additions whose equations are verified in open sources, and two can be fixed only partly. **The best return is in the throat and the boundary layer.** The Kliegel–Levine (KL) re-expansion of Hall's series in ε = 1/(1+R) is about 30 lines added to `hall_uv`. On a numerical check it gives axis Mach 0.790 and Cd 0.9817 at Rc/Rt = 0.625, against measured values of about 0.8 and 0.985. Levine & Coats found it matches data "down to approximately .4", though they call that agreement "fortuitous". **The +25% δ\* bias most likely has one dominant cause**: rockmoc uses a single thickness, while on a cooled wall the thermal layer is much thicker than the velocity layer (φ/θ up to 5–10 at the throat). Adding Bartz's energy-integral equation and two-thickness shape factor, started at the chamber rather than the throat, reduces δ\* by tens of percent, and NASA Lewis measured exactly this mechanism. Relaminarization needs only the published switch K_ax > 2.88×10⁻⁶. Finite-rate chemistry can be built on Cantera from derived quasi-1D equations, but the JANNAF ODK reaction set and start logic were not retrieved. For separation, the criteria (Schmucker, Stark, Morrisette–Goldberg, Kalt–Badal) are well sourced, but internal-shock detection by characteristic crossing is a design inference, not a published algorithm. Two-phase loss is the weakest area: the coefficients of every empirical correlation (Hermsen D43, SPP) are unverified, so only a derived gas–particle ODE is safe to code. Finally, Zebbiche's thesis supplies an open, perfect-gas axisymmetric MLN table. It also shows that the "Sivells Mach 4" output widely used as an axisymmetric benchmark is in fact a planar run.

A note on evidence status, which applies throughout. Every equation below is tagged **VERIFIED**, **RECONSTRUCTED** or **UNVERIFIED**:

- **VERIFIED**: read in a fetched primary or secondary source and, where noted, checked numerically.
- **RECONSTRUCTED**: rebuilt from OCR-damaged scans or derived from conservation laws, and passed internal-consistency checks.
- **UNVERIFIED**: recalled, or seen only as a search snippet. Do not hard-code it until the original has been read.

## KL series fixes the throat down to Rc/Rt ≈ 0.4; below that only a time-marching solver works

Hall's cylindrical expansion in 1/R fails for a structural reason. Levine & Coats carried both cylindrical and true toroidal-coordinate solutions to 5th order and found that "none of them appear to be capable of yielding valid, accurate solutions for small R". Hall's 3rd order **diverges below R ≈ 1.5**, and every higher order diverges below R ≈ 2. Their explanation is that at small R the upstream geometry, not just the local curvature, controls the throat flow ([Levine & Coats 1970, NTRS 19710000970](https://ntrs.nasa.gov/api/citations/19710000970/downloads/19710000970.pdf)). So more terms of any local series will not help. Dutton & Addy's series (TRANNOZ) belongs to the same family of local expansions in 1/R; nothing found shows it works below R = 1 ([Illinois Experts, TRANNOZ](https://experts.illinois.edu/en/publications/trannoz-a-computer-program-for-analysis-of-transonic-throat-flow-/)).

The practical exception is Kliegel & Levine (1969). They introduced no new physics: they re-expanded Hall's own u1–u3 and v1–v3 in ε = 1/(R+1), a parameter that never exceeds 1. Levine & Coats report "excellent agreement with the data for R down to approximately .4". They add that the KL series does not satisfy the cylindrical equations and is "actually an empiricism", but "still … a useful and unique engineering tool". Carrying it to 4th order made it worse, and 5th order diverged ([Levine & Coats 1970](https://ntrs.nasa.gov/api/citations/19710000970/downloads/19710000970.pdf)). The working form, transcribed from a Blida thesis that quotes KL, is ([Kbab 2012, eqs. II.31–II.38](https://di.univ-blida.dz/jspui/bitstream/123456789/3283/1/036-2012.pdf)):

```
ε = 1/(R+1),   z = x·√(2R/(γ+1)),   r = y/y_t
u = 1 + ε u1 + ε²(u1 + u2) + ε³(u1 + 2u2 + u3)
v = √((γ+1)ε/2) · [ ε v1 + ε²(3/2 v1 + v2) + ε³(15/8 v1 + 5/2 v2 + v3) ]

u1 = r²/2 − 1/4 + z                     v1 = r³/4 − r/4 + r z
u2, v2 : identical to rockmoc hall_uv (checked term by term)
u3 = (556γ²+1737γ+3069)/10368 r⁶ − (388γ²+1161γ+1881)/2304 r⁴ + (304γ²+831γ+1242)/1728 r²
     − (2708γ²+7839γ+14211)/82944
     + z[(52γ²+51γ+327)/384 r⁴ − (52γ²+75γ+279)/192 r² + (92γ²+180γ+639)/1152]
     + z²[−(7γ−3)/8 r² + (13γ−27)/144] + (4γ²−57γ+27)/144 z³
v3 = (6836γ²+23031γ+30627)/82944 r⁷ − (3380γ²+11391γ+15291)/13824 r⁵
     + (3424γ²+11271γ+15228)/13824 r³ − (7100γ²+22311γ+30249)/82944 r
     + z[(556γ²+1737γ+3069)/1728 r⁵ − (388γ²+1161γ+1881)/576 r³ + (304γ²+831γ+1242)/864 r]
     + z²[(52γ²+51γ+327)/192 r³ − (52γ²+75γ+279)/192 r] − (7γ−3)/12 r z³
```

This block is **RECONSTRUCTED**. OCR lost the signs, which were restored using rockmoc's existing Hall terms and an irrotationality check. The OCR also printed "1735" in v3's z·r⁵ term; the value 1737 is the one irrotationality requires. With that correction, ∂u/∂r = ∂v/∂x holds to machine precision. The binomial weights are exactly what you get by substituting 1/R = ε/(1−ε) into Hall's series, which confirms how KL was constructed. **One real ambiguity remains.** The OCR'd scaling line could also be read as z = x√(2(R+1)/(γ+1)). Throat-plane values do not depend on the choice, but axial gradients do. The original paper ([AIAA J 7(7):1375, doi 10.2514/3.5355](https://arc.aiaa.org/doi/10.2514/3.5355)) was not accessible, so the scaling must be checked against it or against Cuffel's axial data before release.

The numerical check run for this research (γ = 1.4, throat plane) shows why KL is worth adopting:

| R = Rc/Rt | Method | Axis M | Wall M |
|---|---|---|---|
| 0.625 | Hall 2nd order | 1.300 | 1.114 |
| 0.625 | Hall 3rd order | 0.252 | 2.20 |
| 0.625 | **KL 3rd order** | **0.790** | **1.314** |
| 1.0 | KL | 0.837 | 1.226 |
| 2.0 | KL / Hall 2nd / Hall 3rd | 0.897 / 0.925 / 0.873 | — |

The measured references at R = 0.625 are axis M ≈ 0.8, wall M ≈ 1.4 and Cd = 0.985 (Cuffel, Back & Massier 1969). These figures come from the project brief and were not re-verified against the paywalled paper. KL therefore gets the axis value right and is about 7% low at the wall.

The discharge coefficient has a closed-form series that "converges for all values of R" ([Li et al., FLOMEKO 2019](https://imeko.org/publications/tc9-2019/IMEKO-TC9-2019-028.pdf)). The FLOMEKO source also notes that Hall's own Cd turns negative at small R:

```
Hall (1962):       1 − Cd = (γ+1)/R² · [1/96 − (8γ+21)/(4608R) + (754γ²+1971γ+2007)/(552960R²)]
Hall, KL-corrected: 1 − Cd = (γ+1)/R² · [1/96 − (8γ+21)/(2304R) + (754γ²+2123γ+2553)/(552960R²)]
KL in 1/(1+R):     1 − Cd = (γ+1)/(1+R)² · [1/96 − (8γ−27)/(2304(1+R)) + (754γ²−757γ+3633)/(276480(1+R)²)]
```

Signs were partly reconstructed from the OCR layout. Two cross-checks make these formulas effectively **VERIFIED**. First, Hall's original reproduces the published Cd = 0.99859 at R = 4 ([FLOMEKO 2013](https://imeko.org/publications/tc9-2013/IMEKO-TC9-2013-105.pdf)). Second, the KL formula (0.9817 at R = 0.625) agrees to within 5×10⁻⁴ with an independent integration of the KL velocity field over the throat plane (0.9812). That agreement validates the Cd signs and the u/v transcription together.

KL also has a start-line hazard. At R = 0.625 the KL wall Mach climbs steeply downstream: 1.45 at x/Rt = 0.05, 1.62 at 0.10 and 2.17 at 0.20. A vertical start line at the station where the axis reaches M = 1.02 (x ≈ 0.2) would therefore start the wall too fast. The axis sonic point is predicted at x\*/Rt ≈ +0.173 (R = 0.625), falling to +0.139 (R = 2). The safer choice is a non-vertical initial-value line, running from the axis at M ≈ 1.01 back to the wall near x ≈ 0–0.05, where KL is most trustworthy. This is an inference and has not been tested against data. For R below about 0.4, or for steep 75–90° inlets, the only rigorous option is a numerical transonic solver. Prozan's time-dependent finite-difference method "performed very well" on JPL contours with R = 0.5 and 0.25 ([Prozan, NTRS 19700017856](https://ntrs.nasa.gov/api/citations/19700017856/downloads/19700017856.pdf)). For rockmoc that would be a small numpy MacCormack Euler solver over x/Rt ∈ [−1.5, 0.5], with p0/T0 imposed at the inlet and extrapolation at the supersonic outlet, sized at about 300 lines (estimate).

## The missing thermal thickness explains rockmoc's δ\* bias

Bartz's integral method, the parent of every JPL and NASA Lewis nozzle BL code, carries **two** thicknesses: momentum θ and energy φ. Each has its own 1/7-power profile, over δ and over a separate thermal thickness Δ ([Bartz, NASA CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)). rockmoc already solves the momentum equation (Bartz Eq. 21–22) in the same form. What it lacks is the energy equation (Bartz Eq. 28–30) for an axisymmetric wall:

```
d/dx[r ρ U cp (T0 − Tw) φ] = r q_w
dφ/dx = St·(Taw − Tw)/(T0 − Tw) − φ[(1/ρU) d(ρU)/dx + (1/r) dr/dx − (1/(T0 − Tw)) dTw/dx]
St = (Cf(Re_φ)/2)·(φ/θ)^n / {1 − 5(Cf(Re_φ)/2)^½ [1 − Pr + ln(6/(5Pr+1))]}
```

The structure of these equations is **VERIFIED**. Exactly where the Taw/T0 factor sits in the St term is **OCR-uncertain**. NASA Lewis found that "the best value for the interaction exponent N was determined to be zero" ([Boldman & Graham, TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)). With n = 0 the energy equation decouples and has a closed-form quadrature (Bartz Eq. 44–45). In Bartz's own heated-air comparison, the n = 0 / Coles-adiabatic-Cf variant matched the throat data, while the n = 0.1 / film-temperature variant ran about 30% higher.

Shape factor H = δ\*/θ is then computed by integrating the two-profile model rather than taken from a Crocco closure. Bartz warns that when the thermal and velocity thicknesses are wrongly assumed equal (Δ/δ = 1), this "does have a large effect on δ\*/θ". He also finds that in a cooled throat "φ/θ, and hence Δ/δ, may approach values as high as 5 and, as a consequence, δ\*/θ becomes negative" ([Bartz CR-62615](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)). Back, Massier & Gier measured negative δ\* upstream of the throat ([JPL TR 32-415](https://ntrs.nasa.gov/api/citations/19650010083/downloads/19650010083.pdf)). Lewis BL surveys suggested φ/θ at the throat "possibly as high as 10" ([TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)).

Integrating the Bartz profiles for this research gives the magnitudes. At M = 3 and Tw/T0 = 0.5, H falls from 3.57 (Δ/δ = 1) to 1.40 (Δ/δ = 2) and to −3.2 (Δ/δ = 4). At Tw/T0 = 0.8 it falls from 4.60 to 3.68 to 1.68. The same check showed that rockmoc's current closure, H = H_i·Tw/Te + r(γ−1)/2·M², reproduces the exact single-thickness Crocco–Busemann integral to within 1–2%. **The closure is therefore not wrong; it is solving the wrong problem.** The fix also moves θ in the right direction: a lower H weakens the (H + 2) thinning term in favourable gradients, so θ rises. That addresses rockmoc's 38–43% θ underprediction on the cooled-inlet case. This argument is an inference and has not been tested yet.

The initial condition matters for the same reason. On the Lewis 30°-15° and 60°-15° rigs, changing uncooled inlet length (and so θ at entry) "resulted in a negligible effect on the throat heat transfer". Cooling the inlet, which thickens the thermal layer, produced "pronounced reductions" ([TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)). Graham & Boldman measured directly that "the thermal layer is much thicker than the velocity layer" in the convergent section ([TN D-5356](https://ntrs.nasa.gov/api/citations/19690022491/downloads/19690022491.pdf)). Bartz requires a specified thickness "appropriate to growth in a pipe entrance region of approximately the combustion-chamber length". So rockmoc must start the integration at the chamber/convergent inlet, not at the throat. Initial θ0 and φ0 can be given by the user or estimated from the flat-plate 1/7 law (δ = 0.37x·Re_x^−0.2, θ = 7δ/72; **UNVERIFIED** as quoted here, but textbook-standard), with φ0 > θ0 when the chamber is cooled. Back et al. provide a sensitivity target: increasing the inlet BL from 5% to 25% of inlet radius reduced h by about 10% ([JPL TR 32-415](https://ntrs.nasa.gov/api/citations/19650010083/downloads/19650010083.pdf)).

The Coles skin-friction table Bartz used is **VERIFIED** and directly codable as an interpolation table. It gives (Cf̄·Re_θ̄ → Cf̄) pairs from 2.51 → 0.00590 to 64.8 → 0.00219 ([Bartz CR-62615, Table A-1](https://ntrs.nasa.gov/api/citations/19650013685/downloads/19650013685.pdf)). Its sublayer temperature Ts (Eq. A-3) and the high-Re extrapolation argument could not be read, so rockmoc should keep its Eckert reference-temperature Cf for now. Bartz's own engineering advice is to multiply predicted transonic h by **1.3–1.4, tapering to 1.0 by M ≈ 4**.

Relaminarization has a sourced, one-line switch. Boldman, Schmidt & Gallagher set Re_θ,cr = 360 in an 8°-nozzle momentum relation and obtained ([TN D-4788](https://ntrs.nasa.gov/api/citations/19680024170/downloads/19680024170.pdf)):

```
K_ax = (ν0/u∞²)(du∞/dx) + 0.352·(ν0/(u∞ r))(dr/dx),     K_ax,cr = 0.00368/1278 = 2.88×10⁻⁶   (ν0 at T0)
```

These are **VERIFIED**. Above K_cr, measured throat heat transfer was about **50% of the turbulent prediction**. The Lewis review says that above "about 3×10⁻⁶ a rapid dropoff in the Stanton number occurs", that inlet cooling does not move the onset, and that roughness narrows the laminarizing range ([TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)). The proposed implementation goes beyond what the sources specify. It would switch Cf and St to laminar closures (H_i ≈ 2.2–2.6) when K > K_cr or Re_θ < 360, blend the two over 2.5–4×10⁻⁶, and revert downstream. No validated re-transition criterion exists for the supersonic side. Large high-Pc engines rarely trigger the switch, but the low-p0 Lewis rigs that rockmoc validates against can.

Calibration should be judged against what industry codes achieve. TDK's BL module overpredicts outer-layer thickness and overpredicts heat flux by about 100% at ε = 50, falling to about 3% at ε = 975 (snippet level, [NTRS 19960029260](https://ntrs.nasa.gov/api/citations/19960029260/downloads/19960029260.pdf)). RocketIsp notes that the NASA SP-8120 BL correlation "predicts about twice as much boundary layer loss" as its TDK-trained fit ([RocketIsp](https://rocketisp.readthedocs.io/en/latest/nozzle_eff.html)). rockmoc's present ±20–40% errors are typical of 1960s integral methods. A realistic goal is the ±20% throat heat-transfer accuracy the Lewis energy method achieved ([TN D-6594](https://ntrs.nasa.gov/api/citations/19720008558/downloads/19720008558.pdf)), not CFD-grade accuracy.

## Kinetics is buildable on Cantera; two-phase must stay a derived ODE

Finite-rate chemistry in JANNAF practice comes from ODK (quasi-1D) and TDK, which is "a Two-Dimensional Method of Characteristics solution with fully coupled finite rate kinetics" ([Dunn & Coats, NTRS 19960029260](https://ntrs.nasa.gov/citations/19960029260)). The ODK equation set and its standard reaction table were not retrieved; the NTRS manuals are unreadable scans. The governing system is nonetheless standard and **RECONSTRUCTED (derived)** from conservation laws, with M taken as the *frozen* Mach number:

```
du/dx = u/(1 − M²) · [ −(1/A) dA/dx + Σ_k (W̄/W_k − h_k/(cp T)) dY_k/dx ]
dT/dx = −( u du/dx + Σ_k h_k dY_k/dx ) / cp
dp/dx = −ρ u du/dx,      dY_k/dx = ω̇_k W_k / (ρ u)
```

Every term maps onto Cantera: `net_production_rates`, `partial_molar_enthalpies / molecular_weights`, `cp_mass`, `mean_molecular_weight`. The system is stiff, so integrate it with `solve_ivp(method='BDF'/'LSODA')`. Because ṁ is an eigenvalue at the frozen-sonic point, the cheapest robust start is to use shifting equilibrium (from the existing RocketCEA/Cantera backend) up to A/A\* ≈ 1.01–1.1 on the supersonic side, then integrate kinetically along the MOC contour's A(x). The notes judge this to be the classical ODK approach, but that is unverified. NASA's LSENS solves the same 1-D problem with assigned area or pressure and the LSODE integrator, and is a documented reference to check against ([NTRS 19950008642](https://ntrs.nasa.gov/citations/19950008642)). The output slots into rockmoc's existing frozen/shifting bounds through RocketIsp's definition, fracKin = (Isp_kin − Isp_frozen)/(Isp_eq − Isp_frozen) ([RocketIsp](https://rocketisp.readthedocs.io/en/latest/nozzle_eff.html)).

The size of the loss is worth knowing before investing. RocketIsp's TDK-trained models put kinetic loss at **0.1–10% of Isp**, and "thruster size (Rthrt) is the most important factor" ([RocketIsp](https://rocketisp.readthedocs.io/en/latest/nozzle_eff.html)). For LOX/RP-1 the effect is "very small", while for LOX/LH2 at high area ratio and low MR "CEA gets more than 3 seconds higher than TDK" ([RocketCEA: CEA vs TDK](https://rocketcea.readthedocs.io/en/latest/cea_tdk_rpa.html)). So the kinetic module matters for small, low-Pc, high-ε engines and hydrogen, not for the booster-class designs rockmoc mostly targets.

An intermediate option is Bray's sudden-freezing criterion: integrate in equilibrium until the Lagrangian rate of change of equilibrium composition equals the one-way recombination rate, then freeze. RocketCEA can then expand frozen from that area ratio. Bray's result that flow "will very rapidly 'freeze'" once it departs from equilibrium is **VERIFIED** ([Bray 1959, JFM](https://resolve.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/atomic-recombination-in-a-hypersonic-windtunnel-nozzle/A180D2661B68DD9CB329906B6F33D23D)), as is the caveat that with small expansion angles freezing "stretches out over quite some length" ([Watson 1961, Caltech](https://authors.library.caltech.edu/records/rz0ga-7sp02)). The exact form of the criterion is **UNVERIFIED**; Bray's open-access [Southampton thesis](https://eprints.soton.ac.uk/458375/1/38084.pdf) is the place to confirm it. For coupled multi-reaction rocket exhaust, Bray will bracket ODK rather than match it.

Two-phase loss is the area where the notes find the least verifiable data. The mechanism (velocity and thermal lag of Al2O3 droplets) is well sourced, and the no-lag limit is an upper bound on performance ([Nakka, citing Hoglund 1962](https://www.nakka-rocketry.net/th_2phf.html)). Two items are sourced only secondhand: SPP treats D43 as "a unique function of motor throat diameter", with a rule of thumb of "1μ of droplet diameter … for every inch of nozzle throat diameter" ([US patent 8051640](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8051640)). The Hermsen D43 correlation (3.6304·Dt^0.2932·[1 − exp(−0.0008163·ξ·Pc·τ_res)]) and the SPP loss-correlation exponents are **UNVERIFIED** recollections; AD-A004666 returned 403 and the AIAA originals are paywalled. One 2019 study even found loss *decreasing* with particle size at fixed inert loading ([BIT 2019](https://pure.bit.edu.cn/en/publications/effects-of-particle-size-on-two-phase-flow-loss-in-aluminized-sol)), which shows how much the result depends on regime. The defensible addition is the **RECONSTRUCTED (derived)** gas–particle quasi-1D ODE below. Its no-lag limit should match CEA run with condensed Al2O3:

```
du_p/dx = (u − u_p)/(τ_v u_p) · f_D,     τ_v = ρ_p d²/(18μ),    f_D = C_D Re_p/24
dT_p/dx = (T − T_p)/(τ_T u_p) · Nu/2,    τ_T = ρ_p c_s d²/(12k)
ṁ_g(h + u²/2) + ṁ_p(c_s T_p + u_p²/2) = const
Isp = (ṁ_g u_e + ṁ_p u_pe + p_e A_e)/(ṁ_g + ṁ_p)/g0
```

The particle diameter d stays a user input until Hermsen is verified. RocketIsp itself still asks for a user-supplied η_TP ([RocketIsp](https://rocketisp.readthedocs.io/en/latest/nozzle_eff.html)), so rockmoc would be ahead of that reference tool even without the correlation. For combining losses, RocketIsp documents both a multiplicative stack and a form that subtracts the BL loss, and states that "either form may be correct" depending on how the BL loss was derived ([RocketIsp efficiencies](https://rocketisp.readthedocs.io/en/latest/_sources/efficiencies.rst.txt)). rockmoc computes BL loss physically, which favours Isp = Isp_ODE·η_div·η_kin·η_TP − ΔIsp_BL in `performance/integrated.py`.

## Separation needs a pluggable criterion; internal shocks need a crossing test

Every one-line separation criterion carries large scatter. NASA advised an extra **20% margin** ([Östlund 2004](https://www.mech.kth.se/thesis/2004/phd/phd_2004_jan_ostlund.pdf)), and no source gives an RMS error for any criterion. The **VERIFIED** set, all expressed as p_sep/p_a, is:

| Criterion | Formula | Basis / caveat |
|---|---|---|
| Summerfield | ≈ 0.35–0.4 | Conical, 1950s; called "inadequate" by TP-1207 |
| Kalt & Badal | (2/3)(p0/p_a)^−0.2 | Used by RocketCEA `separated_Cf` |
| Schilling | 0.582 (p0/p_a)^−0.195 contoured; 0.541 (p0/p_a)^−0.136 conical | Stark quotes 0.583 |
| Schmucker | (1.88 M_i − 1)^−0.64 | European industry default |
| Stark (DLR 2005) | π/(3 M_sep) | Fit to enlarged database; no dependence on nozzle type |
| Morrisette–Goldberg | 1.082 − 0.363 M_s + 0.0386 M_s², 2.4 ≤ M_s ≤ 4.5 | Contoured axisymmetric |

Sources: [Östlund 2004](https://www.mech.kth.se/thesis/2004/phd/phd_2004_jan_ostlund.pdf), [Stark AIAA 2005-3940](https://elib.dlr.de/49253/1/AIAA2005-3940.pdf), [NASA TP-1207](https://ntrs.nasa.gov/api/citations/19780022611/downloads/19780022611.pdf), [RocketCEA](https://rocketcea.readthedocs.io/en/stable/ambient_isp.html).

At M = 3 these give 0.40, 0.375, 0.349 and 0.341 respectively (Summerfield, Schmucker, Stark, Morrisette–Goldberg). Stark found that DLR cold-flow data "clearly fall below" the expected trend and that "Schmucker's under predict the separation location" ([Stark 2005](https://elib.dlr.de/49253/1/AIAA2005-3940.pdf)). In practice this means the actual separation sits farther downstream than Schmucker says. rockmoc should therefore report Schmucker and Stark as a band rather than pick one.

TP-1207 adds two **VERIFIED** guards. The nozzle flows full when p_t/p_a ≥ 1 + 0.39·(p_t/p_e). Criteria lose validity in the "end effect" zone A_s/A\* > 0.8(A_e/A\* − 1) + 1 ([TP-1207 Eqs. 8–9](https://ntrs.nasa.gov/api/citations/19780022611/downloads/19780022611.pdf)). The physically split two-step criteria (Zukoski p_i/p_p = (1 + 0.5M_i)^−1, plateau recovery p_p/p_a ≈ 0.9 for FSS) are sourced. Östlund's oblique-shock fit, however, is **RECONSTRUCTED**, and its θ/β labels are probably swapped in the extracted text. Östlund also states that no p_p/p_a model "is currently available for contoured nozzles" ([Östlund 2004](https://www.mech.kth.se/thesis/2004/phd/phd_2004_jan_ostlund.pdf)).

Thrust with separation needs no new physics. Integrate the attached MOC wall pressure up to x_sep and set p_w ≈ p_a downstream, which is equivalent to Cf of the nozzle truncated at ε_sep and evaluated at p_a. RocketCEA's Kalt–Badal example is a ready cross-check: LOX/LH2, Pc = 1000 psia, ε = 40 gives Cf 1.78854 → 1.45502 with ε_sep = 27.32 ([RocketCEA](https://rocketcea.readthedocs.io/en/stable/ambient_isp.html)). The finer model ramps wall pressure from p_p to about p_a on the skirt and adds only a small negative correction.

Restricted-shock separation (RSS) cannot be predicted from any criterion, but its precondition can be detected. "A transition from FSS to RSS can only occur in nozzles featuring an internal shock" ([Östlund 2004](https://www.mech.kth.se/thesis/2004/phd/phd_2004_jan_ostlund.pdf)). In DLR's TICTOP tests, "converging characteristics indicate an internal shock". A TICTOP with only moderate compression still showed RSS in 2 of 16 shutdowns, and RSS side loads were about **3× FSS** ([Frey et al., EUCASS 2017](https://www.eucass.eu/doi/EUCASS2017-359.pdf)). So the rockmoc test (a design inference, not a published algorithm) is a sign change or collapse of the cell Jacobian ∂(x,r)/∂(j,k) between adjacent right-running characteristics. A cheap wall pre-check is dθ_w/dx < 0 together with dp_w/dx > 0. The output is `has_internal_shock`, the origin point, and an "RSS-prone" flag whenever predicted separation lies downstream of the shock. Shock-fitting is out of scope; no open source documents it for TOP contours.

## Zebbiche supplies the axisymmetric MLN table; the popular "Sivells Mach 4" is planar

The Argrow & Emanuel (1988) tables are paywalled and were not found reproduced anywhere ([J Fluids Eng 110(3), doi 10.1115/1.3243546](https://link.springer.com/article/10.1186/s42774-025-00244-1)). The usable open substitute is Zebbiche's Doctorat d'État thesis. Its perfect-gas rows (γ = 1.402, straight sonic line, y\* = throat radius) are **VERIFIED** by transcription and pass a derived A/A\* check: (y_S/y\*)² = 1.6858 at M = 2, against the isentropic 1.6875 ([Zebbiche thesis, Tables 5.42–5.49](https://di.univ-blida.dz/jspui/bitstream/123456789/6224/1/32-620-91-1.pdf)).

| M_S | L/y\* | θ\* corner (°) | θ_max (°) | M\* | y_S/y\* | Inflection (x/y\*, y/y\*) |
|---|---|---|---|---|---|---|
| 1.5 | 2.0143 | 2.4138 | 3.7584 | 1.1515 | 1.0843 | 0.4317, 1.0261 |
| 2.0 | 3.5731 | 5.7466 | 8.2830 | 1.2846 | 1.2984 | 0.5666, 1.0769 |
| 3.0 | 8.3554 | 11.7091 | 15.4652 | 1.4938 | 2.0542 | 0.7074, 1.1845 |
| 4.0 | 16.8144 | 16.1020 | 20.2781 | 1.6427 | 3.2629 | 0.8239, 1.2894 |
| 5.0 | 30.7030 | 19.2722 | 23.5528 | 1.7509 | 4.9748 | — |
| 6.0 | 51.9449 | 21.6141 | 25.8751 | 1.8319 | 7.2440 | — |

Two facts in this table change how rockmoc's tests should be written. First, in an axisymmetric MLN the maximum wall angle is **not** at the corner: the wall keeps turning outward to θ_max at an inflection point near x/y\* ≈ 0.4–0.8. Both θ\* and θ_max are far below ν_e/2. At M = 3, θ\* = 11.7° and θ_max = 15.5° against ν_e/2 = 24.9°. A test that assumes "θ_max at the sharp corner" will fail against this source, and the planar ν_e/2 identity has no axisymmetric counterpart. Second, the thesis's own mesh studies and the 2024 INCAS rerun differ by about 0.5%, so tolerances of **1–2% on L/y\*** are appropriate. Do not use the INCAS Table 1 as a perfect-gas reference. Its A_E/A\* = 1.7711 at M = 2 (perfect gas: 1.6875), and its rows match Zebbiche's high-temperature T0 = 2000 K case ([Abada et al., INCAS 2024](https://bulletin.incas.ro/files/abada_kbab_haif__vol_16_iss_2.pdf)).

For ideal (non-MLN) contours the canonical reference is Sivells AEDC-TR-78-63, whose sample is an axisymmetric Mach 4 nozzle with γ = 1.4, η = 8.67°, RC = 6 and an inviscid exit radius of 12.25 in. The MIT-licensed aldorona/contur repo ships a clean `docs/output.txt` ([aldorona/contur](https://github.com/aldorona/contur)), which evand/nozzle tests quote as "Sivells Mach 4" ([evand/nozzle](https://github.com/evand/nozzle/blob/main/tests/test_contur_equivalence.py)). **That output is almost certainly a planar run.** Its exit-to-throat height ratio, 1.6219617/0.15141471 = 10.712, equals the planar A/A\*(M = 4) = 10.719; an axisymmetric radius ratio would be √10.72 ≈ 3.27. Its SF, BMACH and XD(IN) (7.5526, 3.2, 220.55 in) also differ from the report's axisymmetric sample (21.750386, 3.0821543, 179.7566517 in) ([Sivells.pdf](https://github.com/aldorona/contur/blob/master/sivells/Sivells.pdf)). The file is still a good planar test. To get the axisymmetric case, rerun the F90 port with JD = 0 and BMACH = 3, and confirm that SF = 21.750386 is reproduced. The same report gives one more **VERIFIED** throat-related fact rockmoc can assert: the second derivative of velocity ratio at the sonic point is negative for R < 10.525 in axisymmetric flow.

## Implementation order follows verification, not ambition

The plan below ranks work by verified gain per line of code. Effort figures are estimates.

| # | Change | rockmoc location | Effort | Evidence status | Validation gate |
|---|---|---|---|---|---|
| 1 | `kl_uv` (ε = 1/(1+R) Hall re-expansion) for 0.4 ≤ R < 1.5; KL Cd series for all R; warn below R ≈ 0.4; non-vertical start line | `moc/transonic.py` (keep `hall_uv` for R ≥ 1.5) | ~30–60 lines | Reconstructed, cross-checked (irrotational, Cd ↔ field integral) | R = 0.625: axis M 0.79 ± 0.03, Cd 0.982 ± 0.003 vs 0.985; R = 2 TN D-5377 throat taps unchanged |
| 2 | Zebbiche GP MLN table + test asserting θ\* ≠ θ_max; relabel any planar Sivells test | `validate/experimental.py` | ~50 lines | Verified | L/y\*, θ\*, θ_max, y_S/y\* within 1–2% |
| 3 | Separation criteria (Schmucker, Stark, M–G, Kalt–Badal) + TP-1207 full-flow and end-effect guards + truncated-nozzle Cf | new `analysis/separation.py`, hook into `performance/integrated.py` | ~150 lines | Verified | Reproduce RocketCEA ε_sep = 27.32 / Cf 1.455 case |
| 4 | Bartz energy equation (n = 0), two-thickness H(M, Tw/T0, Δ/δ) cached table, start at chamber inlet with θ0/φ0 | `analysis/boundary_layer.py` | ~200 lines | Structure verified; St factor grouping OCR-uncertain | TN D-5377 δ\* bias < 10%; TN D-6594 Tables V–VII cooled vs adiabatic inlet h(x) within ±20% |
| 5 | K_ax relaminarization switch with laminar Cf/St blend | `analysis/boundary_layer.py` | ~60 lines | Verified threshold; blend is design choice | TN D-4788 low-p0 throat h ≈ 50% of turbulent |
| 6 | Internal-shock detector (characteristic crossing) + RSS-prone flag | MOC marcher | ~80 lines | Design inference | ε = 20 TOP-like contour (r_t 33.54 mm, L 350 mm) flags a shock; TIC contour does not |
| 7 | Quasi-1D finite-rate kinetics (Cantera, equilibrium start at A/A\* ≈ 1.05) → fracKin | new `thermo/kinetics.py` | ~200 lines | Derived; ODK specifics unverified | Lies between frozen and shifting; LOX/RP-1 near-equilibrium; check vs LSENS |
| 8 | Gas–particle 1-D ODE with user d; Hermsen/SPP only after source check | new `performance/two_phase.py` | ~150 lines | Derived; correlations unverified | No-lag limit equals CEA condensed-Al2O3 Isp |
| 9 | 2-D MacCormack transonic Euler for R < 0.4 or steep inlets | `moc/transonic.py` optional | ~300 lines | Method verified (Prozan); code new | Cuffel R = 0.625 wall/axis M; agreement with KL at R = 1–2 |

Items 1–3 rest on verified or cross-checked equations and should come first. Item 2 is the cheapest and protects every later change against a wrong axisymmetric benchmark. Item 4 carries the largest physics payoff but must go through validation on the scanned TN D-6594 tables, which need OCR and manual checking. Before items 7–9 ship, three documents should be retrieved:

- the original KL paper, to settle the z-scaling;
- Bray's thesis, for the exact freezing criterion;
- Hermsen (1981), the SPP manual (AFRPL-TR-83-036) and AD-A004666, for the two-phase coefficients.

Items 8–9 should stay last. Kinetic loss on large LOX/RP-1 engines is "very small", and the time-marching solver only matters below the R ≈ 0.4 range that KL already covers.

## Conclusion

The research changes the diagnosis of two limitations. The δ\* error is not a closure-accuracy problem; rockmoc's Crocco H matches the exact single-layer integral to 1–2%. It is a missing state variable: on cooled walls the thermal layer, not the velocity layer, controls displacement, heat flux and even the sign of δ\*. That makes the energy equation, started in the chamber, the key change in the BL roadmap, while finer shape-factor correlations matter much less. In the same way, the throat limit cannot be fixed by adding terms to the series, which Levine & Coats proved diverges. The fix is a change of expansion parameter that is admittedly empirical, plus an honest numerical fallback below R ≈ 0.4.

The benchmark audit also shows how fragile nozzle validation data is. The open "Sivells Mach 4" axisymmetric reference used by at least one public MOC code appears to be planar, and the planar ν_e/2 corner-angle rule has no axisymmetric equivalent. rockmoc's 30/30 gated checks are only as strong as these references. Pinning the Zebbiche table and a re-run axisymmetric CONTUR case into `validate/experimental.py` before touching the physics is the step that lets every later improvement be measured rather than assumed.
