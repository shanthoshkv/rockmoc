# Validation sources

## `anderson_naca_isentropic`

- **Citation:** Standard isentropic & Prandtl–Meyer tables (γ=1.4); Anderson Fundamentals of Aerodynamics; NACA Report 1135 lineage

- **Type:** analytical_table

- **Notes:** Closed-form perfect-gas relations; reference digits fixed to commonly published table values


## `nasa_tm_x_1502`

- **Citation:** Vanco & Goldman, NASA TM-X-1502 (1968), Computer program for design of two-dimensional supersonic nozzle with sharp-edged throat

- **Type:** nasa_report

- **Local PDF:** `refs/NASA-TM-X-1502.pdf`

- **Notes:** Sample Me=4.05, γ=1.36; θ_max=ν_e/2; y_e=A/A*; Δν accuracy discussion Fig. 3


## `naca_rm_e8j12`

- **Citation:** Shames & Seashore, NACA RM E8J12 (1948), Design data for graphical construction of two-dimensional sharp-edge-throat supersonic nozzles

- **Type:** naca_report

- **Local PDF:** `refs/NACA-RM-E8J12.pdf`

- **Notes:** Graphical construction data M=1.2–10; design rules θ_max=ν_e/2 and y_e/y_t=A/A*


## `nasa_tm_x_2343`

- **Citation:** Goldman & Vanco, NASA TM-X-2343 (1971), Computer program for design of two-dimensional supersonic nozzle with boundary-layer correction

- **Type:** nasa_report

- **Local PDF:** `refs/NASA-TM-X-2343.pdf`

- **Notes:** BL displacement geometry: r_corr = r_inviscid + δ*


## `rao_1958`

- **Citation:** G.V.R. Rao, Exhaust Nozzle Contour for Optimum Thrust, Jet Propulsion / ARS (1958)

- **Type:** journal

- **Local PDF:** `refs/Rao_1958_Exhaust_Nozzle_Optimum_Thrust.pdf`

- **Notes:** TOP length often ~60–90% of 15° conical; θ_n > θ_e engineering correlations


## `lee_vernacchia_spike`

- **Citation:** C.C. Lee (1963) annular aerospike; Vernacchia Spike Contour Algorithm notes

- **Type:** technical_note

- **Local PDF:** `refs/Spike_Contour_Algorithm.pdf`

- **Notes:** δ = π/2 − ν_e; isentropic Prandtl–Meyer plug geometry


## `sutton_1d_cf`

- **Citation:** Sutton & Biblarz, Rocket Propulsion Elements — ideal 1-D thrust coefficient

- **Type:** textbook

- **Notes:** Cf formula with pe/pc, pa/pc, ε, Me perfect gas


## `zucrow_hoffman`

- **Citation:** Zucrow & Hoffman, Gas Dynamics — axisymmetric MOC unit processes

- **Type:** textbook

- **Notes:** C± slopes tan(θ±μ); axisym source −sinθ sinμ /(r cos(θ±μ))

