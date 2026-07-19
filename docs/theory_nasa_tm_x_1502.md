# Theory: NASA TM-X-1502 planar sharp-throat MLN

**Primary source:** Vanco & Goldman, NASA TM-X-1502, NTRS 19680005278 (1968).  
**Local PDF:** `refs/NASA-TM-X-1502.pdf`

## Scope

- Two-dimensional (**planar**) flow only  
- Perfect gas, constant γ, isentropic, inviscid  
- Sharp-edged throat  
- Design goal: **uniform parallel exit** (minimum length)

## Design rules implemented in `design_planar_mln_nasa`

| Rule | Value |
|------|--------|
| Wall angle at throat | \(\theta_{\max} = \nu_e / 2\) |
| Mesh control | \(\Delta\nu\) (degrees input; adjusted to integer intervals) |
| Internal scale | throat half-height \(y_t = 1\) |
| Exit height (uniform exit) | \(y_e / y_t = A/A^*(M_e)\) (1-D isentropic; paper notes exit \(y\) known from 1-D) |
| Characteristic slopes | Family I \(\tan(\mu+\theta)\), Family II \(-\tan(\mu-\theta)\), averaged |

## Accuracy (paper Fig. 3, \(M_e=4.05\), \(\gamma=1.36\))

| \(\Delta\nu\) | Approx. error in exit \(x\) |
|---------------|-----------------------------|
| 0.1° | ~0.5% |
| 0.5° | ~2.5% |

Prefer \(\Delta\nu \le 0.1^\circ\) for production designs. Use `mesh_convergence_planar()`.

## Sample input (paper)

```
ME = 4.05
DV = 0.10   (degrees)
GAM = 1.36
```

```python
from rockmoc import design_planar_mln_nasa
r = design_planar_mln_nasa(4.05, 1.36, dnu_deg=0.1, throat_half_height=0.05)
print(r.summary())
```

## Related NASA work

- **TM-X-2343** (1971): same geometry + boundary-layer correction  
- Edelman Meteor 22; NACA RM E8J12 (tabulated sharp-throat contours)  
- Shapiro, *Dynamics and Thermodynamics of Compressible Fluid Flow*

## Not in this paper

Axisymmetric rocket bells — see axisymmetric MOC designers (separate module). Aerospike — see Lee / Vernacchia.
