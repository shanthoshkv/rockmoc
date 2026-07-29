# Glossary & symbol sheet

Quick reference for rockmoc homework and lab reports.  
API angles are **radians** unless a field name ends in `_deg`.

---

## Greek & Latin symbols

| Symbol | Name | Meaning in rockmoc |
|--------|------|-------------------|
| \(\gamma\) | gamma | Ratio of specific heats \(c_p/c_v\) |
| \(R\) | specific gas constant | J/(kg·K); \(R = R_u / \mathcal{M}\) |
| \(M\) | Mach number | \(V/a\); sonic when \(M=1\) |
| \(M_e\) | exit Mach | Design exit Mach (`exit_mach`) |
| \(A\) | area | Cross-section; throat sonic area \(A^*\) |
| \(\varepsilon\) | area ratio | \(A_e/A^*\); axisym \(\varepsilon=(R_e/R_t)^2\) |
| \(R_t\), \(R_e\) | throat / exit radius | Wall radii [m] |
| \(L\) | nozzle length | Axial span of wall \(x_e - x_t\) |
| \(L/R_t\) | nondimensional length | Primary size metric |
| \(\theta\) | flow / wall angle | Angle from axis [rad] |
| \(\theta_{\max}\) | max wall angle | Peak turning near throat |
| \(\theta_e\) | exit wall angle | Should be ~0 for MLN/ideal design exit |
| \(\nu\) | Prandtl–Meyer angle | Turning from \(M=1\) to local \(M\) |
| \(\nu_e\) | exit PM angle | \(\nu(M_e,\gamma)\) |
| \(\mu\) | Mach angle | \(\mu=\arcsin(1/M)\) |
| \(C^\pm\) | characteristics | Right/left-running Mach waves |
| \(p_c\), \(p_e\), \(p_a\) | pressures | Chamber, exit static, ambient |
| \(p_e/p_c\) | pressure ratio | Isentropic function of \(M_e,\gamma\) |
| \(C_F\) | thrust coefficient | \(F/(p_c A_t)\) |
| \(c^*\) | characteristic velocity | Chamber performance metric |
| \(I_{sp}\) | specific impulse | \(F/(\dot m g_0)\) [s] |
| \(\eta_{\mathrm{div}}\) | divergence efficiency | \(\approx(1+\cos\theta_e)/2\) |
| \(\delta^*\) | displacement thickness | BL wall growth (production) |
| \(\delta\) | geometry flag / spike angle | 1=axisym; aerospike \(\delta=\pi/2-\nu_e\) |

---

## Nozzle types (one line each)

| Type | Code | One-liner |
|------|------|-----------|
| MLN | `design_mln` | Shortest contoured nozzle for given \(M_e\) |
| Ideal | `design_ideal` | Full-length contoured reference |
| TIC | `design_tic` | Truncated ideal (e.g. 80% length) |
| Rao TOP | `design_rao` | Thrust-optimized engineering bell |
| Conical | `design_conical` | Straight wall, textbook baseline (15°) |
| Planar MLN | `design_planar_mln_nasa` | 2D NASA TM-X-1502 sharp-throat |
| Aerospike | `design_spike` | Annular altitude-compensating plug |

Full explainers: `from rockmoc import explainer; print(explainer("mln"))`  
or [nozzle_explainers.md](nozzle_explainers.md).

---

## Common abbreviations

| Abbr | Meaning |
|------|---------|
| MOC | Method of Characteristics |
| MLN | Minimum-Length Nozzle |
| TIC | Truncated Ideal Contour |
| TOP | Thrust-Optimized Parabola (Rao family) |
| TOC | Thrust-Optimized Contour |
| IVL | Initial-Value Line (throat start data) |
| BL | Boundary layer |
| CEA | Chemical Equilibrium with Applications (NASA) |
| CFD | Computational Fluid Dynamics |

---

## Units cheat sheet

| Quantity | rockmoc API | Student report tip |
|----------|-------------|--------------------|
| Length | m | Use `result.report(length_unit="mm")` |
| Angle | rad | Use `angle_unit="deg"` in reports |
| Pressure | Pa | Convert: 1 bar = 1e5 Pa; 1 psi ≈ 6895 Pa |
| Temperature | K | |
| R | J/(kg·K) | |

---

## Important identities (perfect gas)

**Area–Mach**

\[
\frac{A}{A^*}=\frac{1}{M}\left[\frac{2}{\gamma+1}\left(1+\frac{\gamma-1}{2}M^2\right)\right]^{\frac{\gamma+1}{2(\gamma-1)}}
\]

**Prandtl–Meyer**

\[
\nu(M)=\sqrt{\frac{\gamma+1}{\gamma-1}}\arctan\sqrt{\frac{\gamma-1}{\gamma+1}(M^2-1)}-\arctan\sqrt{M^2-1}
\]

**Planar MLN corner**

\[
\theta_{\max}=\nu_e/2
\]

**Cone length (half-angle \(\theta_c\))**

\[
L=\frac{R_e-R_t}{\tan\theta_c},\quad R_e=R_t\sqrt{\varepsilon}
\]

**Ideal vacuum thrust coefficient** (structure)

\[
C_{F,\mathrm{vac}}=\sqrt{\frac{2\gamma^2}{\gamma-1}\left(\frac{2}{\gamma+1}\right)^{\frac{\gamma+1}{\gamma-1}}\left(1-\left(\frac{p_e}{p_c}\right)^{\frac{\gamma-1}{\gamma}}\right)}+\varepsilon\frac{p_e}{p_c}
\]

---

## Code ↔ symbol map

| Code | Symbol |
|------|--------|
| `spec.exit_mach` | \(M_e\) |
| `spec.area_ratio` | \(\varepsilon\) target |
| `spec.throat_radius` | \(R_t\) |
| `spec.gas.gamma` | \(\gamma\) |
| `metrics.area_ratio` | geometric \(\varepsilon\) |
| `metrics.length` | \(L\) |
| `metrics.max_wall_angle` | \(\theta_{\max}\) [rad] |
| `contour.theta[-1]` | \(\theta_e\) [rad] |
| `metrics.Cf_vac_1d` | \(C_{F,\mathrm{vac}}\) |

---

*See also: [how_to_use.md](how_to_use.md), [student_problems.md](student_problems.md).*
