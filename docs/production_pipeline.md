# Production pipeline (A–F)

`design_production()` implements the accuracy stack without CFD:

| Step | Feature | Module |
|------|---------|--------|
| A | Fine mesh + expansion table + BL + altitude | `designers/production.py` |
| B | Strict canceling wall (default, with fallback) | `moc/axisym_hf.py` + `moc/strict_wall.py` |
| C | Field-integrated Cf | `performance/field_cf.py` |
| D | Variable-γ tables + `VariableGammaGas` | `thermo/expansion_table.py`, `variable_gamma.py` |
| E | Rao TOC (ideal-MOC length) + Sauer/Hall throat | `designers/rao_toc.py`, `analysis/throat.py` |
| F | UQ + mesh certificate + expanded validation | `analysis/uq.py`, `validate/cases.py` |

## Not included (by request)

- CFD import of throat IVL

## Usage

```python
from rockmoc import DesignSpec, PerfectGas, design_production, design_rao_toc
from rockmoc import monte_carlo_mln, mesh_certificate, field_integrated_cf

prod = design_production(
    DesignSpec(PerfectGas(1.25), throat_radius=0.05, exit_mach=3.0),
    multi_gamma=(1.25, 1.22, 1.18),
    apply_bl=True,
)
toc = design_rao_toc(DesignSpec(PerfectGas(1.2), throat_radius=0.05, exit_mach=3.5), length_fraction=0.8)
uq = monte_carlo_mln(throat_radius=0.05, exit_mach=3.0, n_samples=20)
```
