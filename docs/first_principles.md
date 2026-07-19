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
Expansion tables integrate local \(d\nu/dM\) with \(\gamma(M)\) (multi-γ or CEA frozen samples).  
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
\frac{d(\theta+\nu)}{dx}=-\frac{\sin\theta\sin\mu}{r\cos(\theta+\mu)}
\quad\text{on C+},\qquad
\frac{d(\theta-\nu)}{dx}=-\frac{\sin\theta\sin\mu}{r\cos(\theta-\mu)}
\quad\text{on C-}
\]
Implemented in `moc/equations.py` with predictor–corrector averages (`moc/mesh.py`).

### Planar MLN (NASA TM-X-1502)
- Sharp throat; wall angle at throat \(\theta_{\max}=\nu_e/2\)
- Mesh controlled by \(\Delta\nu\) (prefer \(\le 0.1^\circ\))
- Uniform parallel exit; \(y_e/y_t=A/A^*(M_e)\) in 2D

### Axisymmetric MLN / ideal / TIC
Kernel + canceling wall; strict unit-process wall preferred (`strict_wall.py`).

## 4. Throat (Sauer / Hall / Kliegel)

First-order axisymmetric small-disturbance field near the geometric throat
(\(\varepsilon=R_t/R_c\)):
\[
M^2-1 \approx \alpha\frac{x}{R_t}+\beta\left(\frac{r}{R_t}\right)^2
\]
with \(\alpha\sim\sqrt{(\gamma+1)\varepsilon}\) and Hall factor for finite \(\varepsilon\).
Seeds MOC on a weakly supersonic IVL — not a full CFD throat.

## 5. Boundary layer (NASA TM-X-2343)

Ideal MOC contour is **grown by displacement thickness** \(\delta^*\):
\[
\mathbf{r}_{\mathrm{final}}=\mathbf{r}_{\mathrm{ideal}}+\delta^*\mathbf{n}
\]
TM-X-2343 uses Cohen–Reshotko (laminar) and Sasman–Cresci (turbulent).
rockmoc uses a compressible integral estimate of \(\delta^*\) with the same geometric rule.

## 6. Rao optimum thrust

Rao (1958): fixed length, ambient pressure, and inlet data → wall maximizing thrust
via control surface on a C+ characteristic (calculus of variations + MOC).

rockmoc provides:
- **TOP** engineering parabola with \(\theta_n,\theta_e\) correlations
- **TOC** length referenced to MOC ideal length (`design_rao_toc`)
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
