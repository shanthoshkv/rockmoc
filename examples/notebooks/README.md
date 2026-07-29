# Student notebooks

| Notebook | Topic |
|----------|--------|
| [01_first_mln.ipynb](01_first_mln.ipynb) | First MLN, report, CSV, optional plot |
| [02_compare_nozzles.ipynb](02_compare_nozzles.ipynb) | MLN / ideal / TIC / Rao / cone table |
| [03_performance_altitude.ipynb](03_performance_altitude.ipynb) | Cf & Isp vs altitude |
| [04_nasa_planar.ipynb](04_nasa_planar.ipynb) | NASA TM-X-1502 planar MLN rules |

```bash
pip install -e ".[plot,dev]"
jupyter notebook examples/notebooks
# or: jupyter lab
```

Regenerate JSON sources if needed:

```bash
py -3.12 examples/notebooks/_build_notebooks.py
```

Also see:

- [docs/how_to_use.md](../../docs/how_to_use.md)
- [docs/student_problems.md](../../docs/student_problems.md)
- [docs/glossary.md](../../docs/glossary.md)
