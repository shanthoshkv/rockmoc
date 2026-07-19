# Theory: Annular aerospike — Lee / Vernacchia

**Primary sources**

1. C.C. Lee, *FORTRAN Programs for Plug Nozzle Design*, MSFC / Brown Engineering, March 1963  
2. M. Vernacchia, *Spike Contour Algorithm* (`refs/Spike_Contour_Algorithm.pdf`)  
3. Open code: [mvernacc/aerospike-nozzle-design-gui](https://github.com/mvernacc/aerospike-nozzle-design-gui)

## Assumptions

- Constant \(\gamma\), constant molar mass  
- Adiabatic **isentropic** expansion  
- Expansion fan from **shroud lip**; spike is the streamline that conserves mass on each Mach-line frustum  
- Cylindrical (annular) geometry  

## Key relations

\[
\nu_e = \nu(M_e),\quad \delta = \frac{\pi}{2} - \nu_e
\]

\[
\frac{h_t}{R_e} = \frac{\varepsilon - \sqrt{\varepsilon(\varepsilon - \sin\delta)}}{\varepsilon\sin\delta}
\]

\[
\phi = \nu_e - \nu(M) + \mu(M)
\]

\[
\frac{R_x}{R_e} = \sqrt{\max\left(0,\, 1 - \frac{1}{\varepsilon}\left[\frac{2}{\gamma+1}\left(1+\frac{\gamma-1}{2}M^2\right)\right]^{\frac{\gamma+1}{2(\gamma-1)}}\sin\phi\right)}
\]

\[
\frac{X_x}{R_e} = \frac{1 - R_x/R_e}{\tan\phi}
\]

**Model limit:** requires \(\nu_e < \pi/2\) (otherwise \(\delta < 0\)). Very large \(\varepsilon\) with low \(\gamma\) is rejected with a clear error.

## Mesh

Default \(N = 10\,000\) Mach samples (same as reference GUI). Unit tests use fewer points for speed.

## Performance

Throat momentum + pressure thrust + \(\int (p-p_a)\,dA\) along the spike (trapezoidal cumulative \(I_{sp}\)).

## Usage

```python
from rockmoc import design_spike, design_spike_truncated

r = design_spike(
    Pc=5e6, Tc=1900, gamma=1.27, molar_mass=20.18,
    Re=0.015, expansion_ratio=12, Pa=101325, n_points=10000,
)
trunc = design_spike_truncated(r, length_fraction=0.2)
```
