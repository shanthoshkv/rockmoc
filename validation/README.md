# rockmoc CFD validation package

Open `rockmoc_CFD_validation_report.pdf` (LaTeX; source in `latex/`) or `index.html` for the full report: method, boundary conditions, meshes, per-case mesh plots, contours, centre-line and wall comparisons, every check, limitations.

- `figures/` per-case and summary plots
- `data/` CSV tables, per-case profiles, raw CFD extracts, mesh polygons, complete OpenFOAM case setups
- `FINDINGS.md` narrative findings and the rockmoc fixes

Rebuild: `cd cfd_validation && py -3.12 make_validation_report.py` (needs `results/raw` from `cfdval.py pull`).
