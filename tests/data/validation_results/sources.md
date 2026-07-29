# Validation sources

## `anderson_compressible`

- **Citation:** Anderson, J. D., Modern Compressible Flow (and standard γ=1.4 isentropic/PM teaching tables)

- **Type:** textbook

- **Notes:** Published digits for ν(M), A/A*, pe/pc at γ=1.4


## `naca_report_1135`

- **Citation:** Ames Research Staff, NACA Report 1135 (1953), Equations, Tables and Charts for Compressible Flow

- **Type:** naca_report

- **Notes:** Canonical compressible-flow tables; same closed-form isentropic/PM relations


## `nasa_tm_x_1502`

- **Citation:** Vanco, M. R. & Goldman, L. J., NASA TM-X-1502 (1968), Computer program for design of two-dimensional supersonic nozzle with sharp-edged throat

- **Type:** nasa_report

- **Local PDF:** `refs/NASA-TM-X-1502.pdf`

- **Notes:** Sample Me=4.05, γ=1.36; θ_max=ν_e/2; y_e=A/A*; Δν mesh guidance Fig. 3


## `naca_rm_e8j12`

- **Citation:** Shames, H. & Seashore, F. L., NACA RM E8J12 (1948), Design data for graphical construction of two-dimensional sharp-edge-throat supersonic nozzles

- **Type:** naca_report

- **Local PDF:** `refs/NACA-RM-E8J12.pdf`

- **Notes:** Graphical construction data; design rules θ_max=ν_e/2, y_e/y_t=A/A*


## `nasa_tm_x_2343`

- **Citation:** Goldman, L. J. & Vanco, M. R., NASA TM-X-2343 (1971), Computer program for design of two-dimensional supersonic nozzle with boundary-layer correction

- **Type:** nasa_report

- **Local PDF:** `refs/NASA-TM-X-2343.pdf`

- **Notes:** BL displacement: r_corr = r_inviscid + δ*


## `rao_1958`

- **Citation:** Rao, G. V. R., Exhaust Nozzle Contour for Optimum Thrust, Jet Propulsion / ARS (1958)

- **Type:** journal

- **Local PDF:** `refs/Rao_1958_Exhaust_Nozzle_Optimum_Thrust.pdf`

- **Notes:** TOP length often ~60–90% of 15° conical; θ_n > θ_e


## `lee_vernacchia_spike`

- **Citation:** Lee, C. C. (1963) annular aerospike; Vernacchia Spike Contour Algorithm notes

- **Type:** technical_note

- **Local PDF:** `refs/Spike_Contour_Algorithm.pdf`

- **Notes:** δ = π/2 − ν_e; isentropic PM plug geometry


## `sutton_biblarz`

- **Citation:** Sutton, G. P. & Biblarz, O., Rocket Propulsion Elements — ideal 1-D thrust coefficient

- **Type:** textbook

- **Notes:** Cf formula with pe/pc, pa/pc, ε, Me


## `zucrow_hoffman`

- **Citation:** Zucrow, M. J. & Hoffman, J. D., Gas Dynamics — axisymmetric MOC unit processes

- **Type:** textbook

- **Notes:** C± slopes tan(θ±μ); axisym source −sinθ sinμ/(r cos(θ±μ))


## `huzel_huang`

- **Citation:** Huzel, D. K. & Huang, D. H., Modern Engineering for Design of Liquid-Propellant Rocket Engines (AIAA), conical nozzle practice

- **Type:** textbook

- **Notes:** 15° half-angle conical baseline; L=(Re−Rt)/tan θ_c; divergence loss concept


## `sauer_kliegel`

- **Citation:** Sauer / Kliegel-type transonic throat initial-value line (engineering IVL)

- **Type:** classic_analysis

- **Local PDF:** `refs/Kliegel_Levine_1969.pdf`

- **Notes:** Slightly supersonic throat start; axis θ≈0; M>1 on IVL


## `shapiro_gas_dynamics`

- **Citation:** Shapiro, A. H., The Dynamics and Thermodynamics of Compressible Fluid Flow — Prandtl–Meyer function form

- **Type:** textbook

- **Notes:** ν(M→∞)=π/2(√((γ+1)/(γ−1))−1); characteristic directions

