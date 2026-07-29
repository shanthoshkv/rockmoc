# Nozzle type mini-explainers

Short “why this shape?” notes for homework and lab reports.

You can also print these in Python:

```python
from rockmoc import explainer
print(explainer("mln"))
print(explainer("conical"))
```

---

### Minimum-Length Nozzle (MLN)

Shortest contoured nozzle for a given design Me (sharp throat + canceling wall).

**Why this shape**

- Sharp throat corner starts a Prandtl–Meyer expansion fan.
- Corner angle θ_max ≈ ν_e/2 (planar rule; mild axisym correction).
- Canceling wall turns the flow until θ_e ≈ 0 at design Me.
- Best when length/mass is critical (upper stages, high ε vacuum).
- Exit is design-uniform only approximately (MOC engineering model).

**Use when:** You want the shortest high-ε bell for a target Mach.

**Compare to:** Shorter than ideal and typical TIC; often shorter than Rao TOP.

---

### Ideal contoured bell

Longer classical contoured nozzle: throat arc + canceling region toward uniform exit.

**Why this shape**

- Circular-arc throat (Rc/Rt) then elongated expansion/canceling region.
- Historically the reference for “perfect” inviscid exit uniformity.
- Longer and heavier than MLN or truncated designs.
- Useful as an upper-bound length reference.

**Use when:** You need a full-length contoured reference or textbook ideal shape.

**Compare to:** Longer than MLN, TIC, and typical Rao TOP at the same Me.

---

### Truncated Ideal Contour (TIC)

Cut an ideal contour to a fraction of its length to save mass.

**Why this shape**

- Start from ideal MOC/contour family, truncate axially (e.g. 80%).
- Exit angle θ_e is usually nonzero → divergence loss vs full ideal.
- Classic length–performance trade used on many flight engines.
- `truncate_fraction` is of ideal length (0.05 < f ≤ 1).

**Use when:** You like the ideal family but cannot afford full length.

**Compare to:** Shorter than ideal; longer or shorter than MLN depending on Me and f.

---

### Rao TOP (thrust-optimized parabola)

Industry engineering contour: high thrust for a fixed length (often ~80% of 15° cone).

**Why this shape**

- Parabolic / TOP correlations approximate Rao’s optimum-thrust idea.
- `length_fraction` is relative to a 15° conical nozzle of the same ε.
- Attachment angle θ_n > exit angle θ_e by design.
- Default choice for many liquid-engine preliminary designs.

**Use when:** You want a compact high-thrust bell used in industry practice.

**Compare to:** Usually shorter than full ideal; not identical to pure MOC MLN.

---

### Conical nozzle

Simplest diverging cone (default half-angle 15°) — textbook baseline.

**Why this shape**

- Straight wall at constant half-angle θ_c.
- Length L = (R_e − R_t) / tan(θ_c) with R_e from 1-D ε(Me).
- Divergence efficiency η_div ≈ (1 + cos θ_c)/2 (not axial exit).
- Use as a baseline before claiming MOC/TOP gains.

**Use when:** Homework baseline or “is MOC worth it?” comparisons.

**Compare to:** Usually longer than Rao TOP at 15° for the same ε; large θ_e.

---

### Planar MLN (NASA TM-X-1502)

2D sharp-throat minimum-length nozzle with uniform exit (paper-faithful).

**Why this shape**

- θ_max = ν_e/2 exactly (planar theory).
- Exit half-height y_e/y_t = A/A*(Me).
- Mesh controlled by Δν in degrees (finer → longer runtime).
- Great for learning MOC and matching published NASA samples.

**Use when:** Course labs on planar MOC or TM-X-1502 reproduction.

**Compare to:** Not axisymmetric — do not compare L/Rt directly to bell nozzles.

---

### Annular aerospike (Lee / Vernacchia)

Altitude-compensating plug nozzle from isentropic spike theory.

**Why this shape**

- Spike contour from Prandtl–Meyer turning; δ = π/2 − ν_e.
- Open wake / ambient pressure affects performance with altitude.
- Truncation + base pressure models for practical plugs.
- Requires ν_e < π/2 (very high ε + high γ can fail).

**Use when:** Exploring altitude compensation or plug geometry.

**Compare to:** Different architecture than closed bells — compare Isp vs altitude, not L/Rt alone.

---

### Production design package

Fine MOC + optional multi-γ + BL δ* + field Cf + altitude sweep.

**Why this shape**

- Wraps MLN/ideal/TIC/Rao with accuracy extras.
- BL correction grows the wall by displacement thickness.
- Field Cf integrates exit-plane momentum/pressure.
- Best “submit this contour” path for projects.

**Use when:** Final lab package or export for CFD.

**Compare to:** Same family as underlying kind, with viscous/performance post-processing.
