# Worked homework-style problems

Use these for practice labs. **Try the problem first**, then check the solution sketch.

Requires:

```bash
pip install -e ".[plot,dev]"
```

Symbols: see [glossary.md](glossary.md).

---

## Problem 1 — Thermodynamics warm-up

**Given:** Perfect gas, \(\gamma=1.4\), \(M_e=3.0\).

**Find:** \(\nu_e\) [deg], \(\varepsilon=A_e/A^*\), \(p_e/p_c\).

### Solution sketch

```python
import math
from rockmoc import PerfectGas

gas = PerfectGas(1.4)
Me = 3.0
nu_e_deg = math.degrees(float(gas.prandtl_meyer(Me)))
eps = float(gas.area_ratio(Me))
pe_pc = float(gas.pe_pc_from_mach(Me))
print(nu_e_deg, eps, pe_pc)
```

**Expected (≈):** \(\nu_e \approx 49.76^\circ\), \(\varepsilon \approx 4.235\), \(p_e/p_c \approx 0.0272\).

**Checkpoint:** Errors vs Anderson/NACA tables should be ≪ 0.05% (see deep validation).

---

## Problem 2 — Planar MLN design rule

**Given:** Planar MLN, \(\gamma=1.4\), \(M_e=3\).

**Show:** \(\theta_{\max}=\nu_e/2\), and exit half-height ratio equals \(\varepsilon\).

### Solution sketch

```python
import math
from rockmoc import PerfectGas, design_planar_mln_nasa

gas = PerfectGas(1.4)
Me = 3.0
nu_e = float(gas.prandtl_meyer(Me))
r = design_planar_mln_nasa(Me, 1.4, dnu_deg=1.0, throat_half_height=1.0)
print("θ_max deg", math.degrees(r.theta_max_rad))
print("ν_e/2 deg", math.degrees(0.5 * nu_e))
print("y_e", r.contour.r[-1], "A/A*", float(gas.area_ratio(Me)))
```

**Expected:** \(\theta_{\max}\) matches \(\nu_e/2\) to high precision; \(y_e=\varepsilon\) for \(y_t=1\).

---

## Problem 3 — Axisymmetric MLN package

**Given:** \(R_t=50\) mm, \(\gamma=1.25\), \(M_e=3\), modest mesh.

**Tasks:**

1. Design an MLN with `DesignSpec.student`.
2. Report \(\varepsilon\), \(L/R_t\), \(\theta_{\max}\), \(\theta_e\) in mm/deg.
3. Write a lab-report skeleton.

### Solution sketch

```python
from rockmoc import DesignSpec, design_mln

spec = DesignSpec.student(Me=3.0, gamma=1.25, Rt_mm=50, n=40)
mln = design_mln(spec, store_net=True)
print(mln.report(length_unit="mm", angle_unit="deg"))
mln.lab_report("lab_p3_mln.md", student_name="Your Name", course="AE-4XX")
mln.export_csv("lab_p3_mln.csv")
```

**Discussion prompts:**

- Is \(\theta_e\) near zero? Why does that matter for \(\eta_{\mathrm{div}}\)?
- What happens to \(L/R_t\) if you raise `n` from 40 to 80?

---

## Problem 4 — Cone vs MLN vs Rao

**Given:** Same as Problem 3.

**Tasks:**

1. Design MLN, Rao TOP (80% of 15° cone), and 15° conical.
2. Build a Markdown comparison table.
3. Which is shortest? Which has the largest \(\theta_e\)?

### Solution sketch

```python
from rockmoc import (
    DesignSpec, design_mln, design_rao, design_conical,
    compare_nozzle_designs,
)

spec = DesignSpec.student(Me=3.0, gamma=1.25, Rt_mm=50, n=32)
cmp = compare_nozzle_designs(spec, rao_length_fraction=0.80, include_conical=True)
print(cmp.summary_table())
print(cmp.to_markdown())
cmp.to_csv("lab_p4_compare.csv")
```

**Expected trends:**

- **Conical 15°:** large \(\theta_e\), solid baseline length.
- **Rao 80%:** length near \(0.8\times L_{15^\circ\text{ cone}}\).
- **MLN:** short contoured wall, \(\theta_e\approx 0\).

---

## Problem 5 — Effect of \(\gamma\) at fixed \(M_e\)

**Given:** \(M_e=3\), \(R_t=0.05\) m; compare \(\gamma=1.2\) vs \(1.4\).

**Find:** How \(\varepsilon\) and MLN \(L/R_t\) change.

### Solution sketch

```python
from rockmoc import PerfectGas, DesignSpec, design_mln

for g in (1.2, 1.4):
    gas = PerfectGas(g)
    spec = DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=40)
    r = design_mln(spec)
    print(g, "eps", r.metrics.area_ratio, "L/Rt", r.metrics.length / 0.05)
```

**Expected:** At fixed Me, **higher γ → larger ε** (and usually a longer nozzle).

---

## Problem 6 — TIC trade

**Given:** Ideal and TIC at 80% length, \(\gamma=1.25\), \(M_e=3\), \(R_t=0.05\).

**Show:** \(L_{\mathrm{TIC}} < L_{\mathrm{ideal}}\), and \(\theta_e\) grows when truncated.

### Solution sketch

```python
from rockmoc import PerfectGas, DesignSpec, design_ideal, design_tic
import math

gas = PerfectGas(1.25)
base = dict(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=40)
ideal = design_ideal(DesignSpec(**base))
tic = design_tic(DesignSpec(**base, truncate_fraction=0.8))
print("L", ideal.metrics.length, tic.metrics.length)
print("θ_e deg", math.degrees(abs(ideal.contour.theta[-1])),
      math.degrees(abs(tic.contour.theta[-1])))
```

---

## Problem 7 — Example gases + limits

**Tasks:**

1. List teaching gas presets.
2. Design MLN with `example_gas("lox_rp1")`.
3. Quote known limitations in your report.

### Solution sketch

```python
from rockmoc import example_gas, describe_example_gases, DesignSpec, design_mln, KNOWN_LIMITS

print(describe_example_gases())
gas = example_gas("lox_rp1")
spec = DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.2, n_characteristics=40)
r = design_mln(spec)
print(r.report())
print("Limits:")
for line in KNOWN_LIMITS:
    print(" -", line)
```

---

## Problem 8 — Performance vs altitude (optional)

**Given:** An MLN at \(p_c=5\) MPa, \(T_c=3200\) K.

**Plot / tabulate** \(C_F\) and \(I_{sp}\) from sea level to 80 km.

### Solution sketch

```python
from rockmoc import DesignSpec, design_mln, example_gas, off_design_altitude_sweep

gas = example_gas("lox_rp1")
spec = DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=40)
mln = design_mln(spec)
rows = off_design_altitude_sweep(
    mln.contour, gas, throat_radius=0.05, exit_mach=3.0, pc=5e6, Tc=3200.0
)
for row in rows:
    print(f"{row['altitude_m']/1000:5.0f} km  Cf={row['Cf']:.4f}  Isp={row['Isp_s']:.1f} s")
```

**Discussion:** Why does \(C_F\) rise with altitude for a fixed nozzle?

---

## Grading-friendly checklist (any problem)

- [ ] Units stated (m or mm; degrees when quoting angles)
- [ ] Exactly one expansion target in `DesignSpec`
- [ ] Table or `result.report()` metrics
- [ ] At least one comparison (cone / other type / other γ)
- [ ] Limitations stated (inviscid MOC, frozen γ)
- [ ] CSV or figure attached for the wall

Generate a full Markdown skeleton anytime:

```python
result.lab_report("my_lab.md", student_name="…", course="…")
```
