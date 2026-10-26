# rockmoc Technical Report (LaTeX)

All-in-one reference: purpose, install, CLI/TOML, theory, every public function, examples,
every project figure, validation, limits, changelog.

| File | Description |
|------|-------------|
| [`rockmoc_technical_report.tex`](rockmoc_technical_report.tex) | Source (self-contained bibliography) |
| [`rockmoc_technical_report.pdf`](rockmoc_technical_report.pdf) | Built report (15 pp) |

Outline: 1 Overview · 2 Theory (thermo, multi-γ isentrope, MOC, Hall/Sauer throat, families,
integral BL, loss stack) · 3 API reference · 4 Examples · 5 Figures (all 52 PNGs) ·
6 Validation · 7 Limits · 8 Changelog · References.

Aligned with **0.7.0 + v0.8 accuracy upgrades (unreleased)**.

## Build

Figures load from `../../website/{img,assets,linkedin}`. Regenerate them first if needed:

```powershell
py -3.12 website\make_images.py; py -3.12 website\make_images_extra.py; py -3.12 website\make_assets.py
cd docs\report
pdflatex rockmoc_technical_report.tex; pdflatex rockmoc_technical_report.tex
```

MiKTeX / TeX Live with `amsmath`, `booktabs`, `longtable`, `subcaption`, `listings`, `natbib`,
`hyperref`. Not a flight-qualification document.
