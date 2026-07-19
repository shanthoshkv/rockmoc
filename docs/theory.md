# Theory notes — rockmoc

## Method of characteristics (axisymmetric)

For steady, inviscid, irrotational, isentropic **supersonic** flow, the Euler equations are hyperbolic. Information propagates along **Mach lines** (characteristics):

\[
\left(\frac{\mathrm{d}r}{\mathrm{d}x}\right)_{\pm} = \tan(\theta \pm \mu), \quad \mu = \arcsin\frac{1}{M}.
\]

Along those curves, compatibility relations link the flow angle \(\theta\) and the Prandtl–Meyer function \(\nu(M)\). In **planar** flow the Riemann invariants \(\theta+\nu\) and \(\theta-\nu\) are constant on C+ and C− respectively. In **axisymmetric** flow, geometric source terms \(\propto \sin\theta\sin\mu / r\) appear.

`rockmoc` integrates the unit processes (interior, axis, wall) with predictor–corrector averaging (Zucrow & Hoffman style).

## Prandtl–Meyer function

\[
\nu(M) = \sqrt{\frac{\gamma+1}{\gamma-1}}\arctan\sqrt{\frac{\gamma-1}{\gamma+1}(M^2-1)} - \arctan\sqrt{M^2-1}.
\]

## Contour families

### Minimum-length nozzle (MLN)

Sharp corner at the throat launches a centered expansion fan. The wall then turns back so that expansion waves are canceled and the exit flow is (approximately) axial. For **planar** flow, \(\theta_{\max} = \nu_e/2\). Axisymmetric designs iterate \(\theta_{\max}\) so the kernel axis Mach matches the design \(M_e\).

### Ideal contour

Circular-arc throat, full canceling region sized so the exit core is uniform and parallel at \(M_e\). Longer than MLN for the same expansion ratio.

### Truncated ideal (TIC)

Ideal contour cut at a fraction of ideal length — industry-relevant length/mass trade at some loss of exit uniformity and \(C_F\).

## Performance

Reported \(C_F\) in v0.1 is the **ideal 1D isentropic** thrust coefficient at the design \(M_e\) and geometric area ratio. MOC-integrated exit-plane thrust can be added later from the full field.
