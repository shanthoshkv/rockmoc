# First principles (rockmoc)

This document ties governing equations to code. Local PDFs live in `refs/`.

## 1. Flow assumptions

Steady, inviscid (except post BL correction), adiabatic, **isentropic**, irrotational
supersonic flow of a perfect gas (or table-based isentropic expansion). Gravity neglected.

## 2. Thermodynamics

### Perfect gas
\[
\frac{T_0}{T}=1+\frac{\gamma-1}{2}M^2,\quad
\frac{p_0}{p}=\left(\frac{T_0}{T}\right)^{\gamma/(\gamma-1)},\quad
\frac{A}{A^*}=\frac{1}{M}\left[\frac{2}{\gamma+1}\left(1+\frac{\gamma-1}{2}M^2\right)\right]^{\frac{\gamma+1}{2(\gamma-1)}}
\]

### Prandtl–Meyer
\[
\nu(M)=\sqrt{\frac{\gamma+1}{\gamma-1}}\arctan\sqrt{\frac{\gamma-1}{\gamma+1}(M^2-1)}-\arctan\sqrt{M^2-1}
\]
and \(d\nu=\sqrt{M^2-1}\,dV/V\) along an isentrope.

### Variable γ (along the nozzle)
`build_multi_gamma_table` integrates the thermally perfect isentrope from stagnation
(\(dh=c_p dT\), \(d\ln p=c_p dT/RT\), \(V^2=2\Delta h\), \(A/A^*=(\rho V)^*/\rho V\),
\(d\nu=\sqrt{M^2-1}\,dV/V\)) with \(\gamma(M)\) blended \(\gamma_c\to\gamma_t\to\gamma_e\);
constant γ reproduces perfect gas to 1e-7. CEA frozen samples via `build_cea_frozen_table`.  
`VariableGammaGas` feeds the **same** PM inverse into every MOC unit process, so  
\(\gamma=\gamma(M)\) varies along characteristics and along the wall.  

```python
from rockmoc.thermo import build_multi_gamma_table, VariableGammaGas
from rockmoc import DesignSpec, design_mln
tab = build_multi_gamma_table(gamma_c=1.28, gamma_t=1.24, gamma_e=1.18, R=360.0)
vg = VariableGammaGas(tab)
r = design_mln(DesignSpec(gas=vg, throat_radius=0.05, exit_mach=3.0))
gamma_x = vg.gamma_along_length(r.contour.x)  # γ decreases with x
```

See `thermo/expansion_table.py`, `variable_gamma.py`, `moc/mesh.py` (`_sync_mach`).

## 3. Method of characteristics

### Directions
\[
\left(\frac{dr}{dx}\right)_{\pm}=\tan(\theta\pm\mu),\quad \mu=\arcsin\frac{1}{M}
\]

### Planar (δ=0)
\[
\theta+\nu=\mathrm{const\ on\ C+},\qquad \theta-\nu=\mathrm{const\ on\ C-}
\]

### Axisymmetric (δ=1)
\[
\frac{d(\theta-\nu)}{dx}=-\frac{\sin\theta\sin\mu}{r\cos(\theta+\mu)}
\quad\text{on C+},\qquad
\frac{d(\theta+\nu)}{dx}=+\frac{\sin\theta\sin\mu}{r\cos(\theta-\mu)}
\quad\text{on C-}
\]
(signs verified against exact conical source flow, v0.7).
Implemented in `moc/equations.py` with predictor–corrector averages (`moc/mesh.py`).

### Planar MLN (NASA TM-X-1502)
- Sharp throat; wall angle at throat \(\theta_{\max}=\nu_e/2\)
- Mesh controlled by \(\Delta\nu\) (prefer \(\le 0.1^\circ\))
- Uniform parallel exit; \(y_e/y_t=A/A^*(M_e)\) in 2D

### Axisymmetric MLN / ideal / TIC
Kernel + canceling wall; strict unit-process wall preferred (`strict_wall.py`).

## 4. Throat (Hall / Sauer)

Small-disturbance series in \(\varepsilon=R_t/R_c\) on the vertical line where the axis
reaches `throat_mach`:
- **Hall (1962) 2nd order** — default for axisymmetric \(R_c/R_t\ge1.5\) (`moc/transonic.hall_line`).
  Residual \(\propto\varepsilon^{3.4}\); wall \(p/p_0\) at \(R_c/R_t=2\) within 0.3 % of NASA TN D-5377.
  Diverges for tight throats (axis M 1.30 at 0.625 vs 0.8 measured, Cuffel 1969).
- **Kliegel–Levine (1969)** — Hall's 3rd order re-expanded in \(1/(1+R)\), for \(0.4\le R_c/R_t<1.5\)
  (`kl_line`, `discharge_coefficient`); Cuffel R = 0.625: axis M 0.79 vs 0.8, Cd 0.982 vs 0.985.
- **2-D Euler throat solver** (`moc/transonic_euler`) below 0.4 or for steep inlets.
- **Sauer (1947) 1st order** — planar, and fallback (\(\propto\varepsilon^{2.5}\)).
- Uniform near-sonic line — last fallback (warned).
Kliegel–Levine not implemented. Not a CFD throat.

## 5. Boundary layer (momentum integral, NASA TM-X-2343 geometry rule)

`correct_wall_bl(method="integral")` (default) solves
\[
\frac{d\theta}{ds}=\frac{C_f}{2}-\theta\Big[\frac{2+H-M^2}{1+\frac{\gamma-1}{2}M^2}\frac{1}{M}\frac{dM}{ds}+\frac{\delta}{r}\frac{dr}{ds}\Big]
\]
with Eckert/Sommer–Short reference temperature \(C_f\) (0.0256 Re_θ^{-1/4} turbulent,
0.441/Re_θ laminar), Crocco–Busemann \(H=H_iT_w/T_e+r\frac{\gamma-1}{2}M^2\), Bartz
\(\mu\propto\mathcal M^{0.5}T^{0.6}\), \(Pr=4\gamma/(9\gamma-5)\), adiabatic or prescribed
`wall_temperature`. If the ODE fails → flat-plate closed form of the same closure
(`method="empirical"`, warned); both agree at zero pressure gradient. \(C_f\) within 3 % of
Van Driest II (M ≤ 5). The wall is grown by \(\delta^*=H\theta\):
\[
\mathbf{r}_{\mathrm{final}}=\mathbf{r}_{\mathrm{ideal}}+\delta^*\mathbf{n}
\]
Friction loss \(\Delta C_F=\int\tau_w\cos\alpha\,dA/(p_cA_t)\) → `bl.delta_cf_friction`,
subtracted in `integrated_performance` / `field_integrated_cf` / `design_production`.

## 6. Rao optimum thrust

Rao (1958): fixed length, ambient pressure, and inlet data → wall maximizing thrust
via control surface on a C+ characteristic (calculus of variations + MOC).

rockmoc provides:
- **TOP** engineering parabola with \(\theta_n,\theta_e\) correlations
- **TOC-length TOP** (`design_rao_toc`): TOP parabola scaled to the MOC ideal length
- **TOC by MOC control surface** (`design_rao_control`); \((\theta,M)\) on DE from a 2×2 Newton
  (the θ = μ crossing near L/L15 ≈ 0.8 is regular). Matches the Rao-code bell of NASA
  TM-107285 within 0.37 % of \(r_e\) and SP-8120 chart angles within ~1°.
- **Control-surface TOC** (`design_rao_control`): Rao exit-angle condition  
  \(\tan\theta_e = \sin\mu\cos\mu\big/\big[(\gamma M^2)/(\gamma M^2-1)-\sin^2\mu\big]\)  
  plus MOC kernel + circular arc / Hermite wall (Rao 1958 engineering form)

## 7. Aerospike (Lee / Vernacchia)

Centered expansion from shroud lip; mass continuity on Mach-line frustums fixes
\(R_x(M)\). Truncation + base-pressure and primary cowl are engineering add-ons.

## 8. Performance

- 1-D isentropic \(C_F\)
- Divergence efficiency from exit wall angle
- **Field integral** on exit ray: \(\int(\rho u^2+p)\,2\pi r\,dr\) when MOC net exists
