# rockmoc Technical Report (LaTeX)

Industry-style technical report covering the entire **rockmoc** project: governing equations, contour designers, performance, architecture, multi-source validation, and bibliography.

| File | Description |
|------|-------------|
| [`rockmoc_technical_report.tex`](rockmoc_technical_report.tex) | Full report source (self-contained; no external `.bib`) |

## Contents (outline)

1. Introduction & scope  
2. Flow assumptions  
3. Thermodynamics (isentropic, Prandtl–Meyer, \(C_F\), multi-γ)  
4. Method of characteristics (planar & axisymmetric)  
5. Throat IVL (Sauer / Kliegel)  
6. Designers (MLN, ideal, TIC, Rao, conical, aerospike)  
7. Boundary-layer \(\delta^*\) correction  
8. Performance metrics & field \(C_F\)  
9. Software architecture & API  
10. Validation (43 cases, 12 sources)  
11. Limits & recommended use  
12. References + appendices (module map, equation sheet)

Aligned with package version **0.6.1**.

## Build PDF

Requires a full TeX distribution (TeX Live, MiKTeX, or MacTeX) with common packages: `amsmath`, `booktabs`, `hyperref`, `siunitx`, `natbib`, `listings`, etc.

```bash
cd docs/report
pdflatex rockmoc_technical_report.tex
pdflatex rockmoc_technical_report.tex
```

Or:

```bash
latexmk -pdf rockmoc_technical_report.tex
```

Output: `rockmoc_technical_report.pdf` in this folder.

### Windows (MiKTeX / TeX Live)

```powershell
cd docs\report
pdflatex rockmoc_technical_report.tex
pdflatex rockmoc_technical_report.tex
```

If a package is missing, MiKTeX will usually prompt to install it on first build.

## Notes

- Bibliography is embedded (`thebibliography`); no BibTeX run required.  
- Primary literature PDFs for offline audit live in repo `refs/`.  
- Validation numbers quoted in the report match the stored campaign under `tests/data/validation_results/`.  
- This document is **not** a flight-qualification report; it documents the library for education and preliminary design.
