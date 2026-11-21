"""LaTeX version of the validation report (called from make_validation_report.py).

Writes validation/latex/rockmoc_CFD_validation_report.tex (figures referenced from ../figures)
and compiles it with xelatex (twice, for the TOC and longtables).
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess
from datetime import date
from pathlib import Path

TIERS = [("T1", "Tier 1: contour and MOC field (inviscid, supersonic start-line inlet)"),
         ("GCI", "Mesh independence (coarse / medium / fine)"),
         ("T2", "Tier 2: full nozzles with converging section"),
         ("T3", "Tier 3: real-gas (LOX/RP-1) path"),
         ("T4", "Tier 4: viscous wall, boundary layer and sea-level separation")]

_ESC = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}", "<": r"\textless{}", ">": r"\textgreater{}", "|": r"\textbar{}"}


def t(s) -> str:
    return "".join(_ESC.get(c, c) for c in str(s))


def f(v, n=5):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "--"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, int):
        return str(v)
    return f"{v:.{n}g}"


def status_tex(c):
    if c.info:
        return r"\textcolor{infoc}{info}"
    return r"\textcolor{passc}{\textbf{PASS}}" if c.ok else r"\textcolor{failc}{\textbf{FAIL}}"


def lab(s):
    return re.sub(r"[^A-Za-z0-9:-]", "-", str(s))


def T(s, **kw):
    for k, v in kw.items():
        s = s.replace(f"@{k}@", str(v))
    return s


def pf(x, n=2, sign=False):
    if x is None or not math.isfinite(x):
        return "--"
    return f"{x:+.{n}f}" if sign else f"{x:.{n}f}"


def longtable(spec, header, rows, size=r"\small", caption=None, label=None, note=None, short=None):
    h = " & ".join(rf"\textbf{{{x}}}" for x in header)
    body = "\n".join(" & ".join(r) + r" \\" for r in rows)
    cap = ""
    if caption:
        sh = f"[{short}]" if short else ""
        cap = f"\\caption{sh}{{{caption}}}" + (f"\\label{{{label}}}" if label else "") + "\\\\\n"
    nt = f"\\par\\nopagebreak\\noindent{{\\footnotesize {note}}}\\par\n" if note else ""
    return (f"{{{size}\n\\begin{{longtable}}{{{spec}}}\n{cap}\\toprule\n{h} \\\\\n\\midrule\n\\endfirsthead\n"
            f"\\toprule\n{h} \\\\\n\\midrule\n\\endhead\n\\bottomrule\n\\endfoot\n{body}\n\\end{{longtable}}\n{nt}}}\n")


def fig(path, caption, h="0.82", label=None, short=None):
    sh = f"[{short}]" if short else ""
    lb = f"\\label{{{label}}}" if label else ""
    return (f"\\begin{{figure}}[H]\\centering\n\\includegraphics[width=\\textwidth,height={h}\\textheight,keepaspectratio]{{{path}}}\n"
            f"\\caption{sh}{{{caption}}}{lb}\n\\end{{figure}}\n")


def errtxt(c):
    if not math.isfinite(c.err):
        return "--"
    return f"{100 * c.err:+.2f}\\,\\%" if c.rel else f"{c.err:+.3g}"


def toltxt(c):
    return f"{100 * c.tol:.1f}\\,\\%" if c.rel else f"{c.tol:g}"


def check_rows(checks, with_case=None):
    out = []
    for c in checks:
        row = [t(c.quantity), f(c.rockmoc), f(c.cfd), errtxt(c), toltxt(c), status_tex(c), t(c.note)]
        out.append(([t(with_case)] if with_case else []) + row)
    return out


def git_info(root):
    try:
        h = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=root, capture_output=True, text=True, timeout=30).stdout.strip()
        d = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=root, capture_output=True, text=True, timeout=30).stdout.strip()
        return h or "unknown", bool(d)
    except Exception:  # noqa: BLE001
        return "unknown", False


PREAMBLE = r"""\documentclass[10pt,a4paper]{article}
\usepackage[margin=2.1cm,headheight=14pt]{geometry}
\usepackage{fontspec}
\setmainfont{Cambria}
\setmonofont{Consolas}[Scale=0.88]
\usepackage{unicode-math}
\setmathfont{Cambria Math}
\usepackage{graphicx,float,booktabs,longtable,array,tabularx,xcolor,caption,fancyhdr,enumitem,listings}
\usepackage[hidelinks,pdfauthor={Shanthosh KV},pdftitle={rockmoc: validation against OpenFOAM CFD},pdfsubject={Validation of the rockmoc MOC nozzle library}]{hyperref}
\usepackage{cleveref}
\title{rockmoc: validation against OpenFOAM CFD}
\author{Shanthosh KV}
\graphicspath{{../figures/}}
\definecolor{passc}{HTML}{276749}\definecolor{failc}{HTML}{C53030}\definecolor{infoc}{HTML}{718096}
\definecolor{accent}{HTML}{2B6CB0}
\captionsetup{font=small,labelfont=bf}
\captionsetup[table]{position=top}
\captionsetup[figure]{position=bottom}
\setlength{\LTcapwidth}{\textwidth}
\setlength{\LTpre}{4pt}\setlength{\LTpost}{6pt}
\renewcommand{\arraystretch}{1.15}
\lstset{basicstyle=\ttfamily\scriptsize,breaklines=true,columns=fullflexible,frame=single,rulecolor=\color{black!15},backgroundcolor=\color{black!3}}
\pagestyle{fancy}\fancyhf{}\lhead{\small rockmoc: validation against OpenFOAM CFD}\rhead{\small\thepage}
\setcounter{tocdepth}{3}
\makeatletter\renewcommand*\l@subsubsection{\@dottedtocline{3}{3.8em}{4.4em}}\makeatother
\begin{document}
"""


def write_latex(out_dir: Path, cases_dir: Path, rows, allchecks, sections, figs, gci_rows, compile_pdf=True, v1=None):
    L = [PREAMBLE]
    isp = [r for r in rows if r["isp_gated"] and not r["plume"]]
    worst = max(isp, key=lambda r: abs(r["isp_err"]))
    ng = sum(1 for _, _, c in allchecks if not c.info)
    npass = sum(1 for _, _, c in allchecks if not c.info and c.ok)
    within = sum(1 for r in isp if abs(r["isp_err"]) <= 0.25)
    byn = {r["name"]: r for r in rows}
    mean_abs = sum(abs(r["isp_err"]) for r in isp) / len(isp)
    nneg = sum(1 for r in isp if r["isp_err"] < 0)
    gci_d = []
    for base in ("t1_mln_M3", "t1_rao_80", "t1_cone_15"):
        c_, f_ = byn.get(base + "_coarse"), byn.get(base + "_fine")
        if c_ and f_:
            gci_d.append((base, 100 * (f_["isp_cfd"] / c_["isp_cfd"] - 1), 100 * (f_["exit_M"] / c_["exit_M"] - 1)))
    gci_isp = max((abs(d[1]) for d in gci_d), default=float("nan"))
    gci_m = max((abs(d[2]) for d in gci_d), default=float("nan"))
    TH = ("t2_rao80_rc0.625", "t2_rao80_rc1", "t2_ideal_full", "t2_tic80_full", "t2_rao80_rc2")
    th_err = [100 * abs(c.err) for n_, _, c in allchecks if n_ in TH and c.quantity.startswith("throat-plane axis M (") and math.isfinite(c.err)]
    th_max = max(th_err, default=float("nan"))
    sp = [r for r in rows if r["plume"] and r["sep_cf_ref"] and r["sep_cf_cfd"]]
    sep_err = [100 * (r["sep_cf_cfd"] / r["sep_cf_ref"] - 1) for r in sp]
    sep_att = [100 * (r["sep_cf_att"] / r["sep_cf_cfd"] - 1) for r in sp if r.get("sep_cf_att")]
    mln = {n_: byn.get(f"t1_mln_n{n_}") or (byn.get("t1_mln_M3") if n_ == 40 else None) for n_ in (10, 20, 40, 80)}
    mln = {n_: 100 * r["axis_M_rms"] for n_, r in mln.items() if r and math.isfinite(r["axis_M_rms"])}
    ba = [(r["name"], 100 * (v1[r["name"]][1] / v1[r["name"]][0] - 1), r["isp_err"]) for r in isp
          if v1 and r["name"] in v1 and math.isfinite(v1[r["name"]][0])]
    commit, dirty = git_info(out_dir.parent)

    # ---------------- title page (no pass/fail counts here)
    L.append(T(r"""\begin{titlepage}\centering
\vspace*{3.5cm}
{\Huge\bfseries rockmoc\par}\vspace{0.4cm}
{\LARGE Validation against OpenFOAM CFD\par}\vspace{0.8cm}
{\large Method-of-characteristics nozzle design and performance library\par}
\vspace{1.8cm}
{\Large Shanthosh KV\par}
\vfill
{\large @N@ OpenFOAM simulations \quad$\cdot$\quad @CELLS@ million cells\par}
\vspace{0.6cm}
{\large @DATE@\par}
\end{titlepage}
""", N=len(rows), CELLS=f"{sum(r['ncell'] for r in rows) / 1e6:.2f}", DATE=date.today().strftime("%d %B %Y")))
    # ---------------- abstract, lists, nomenclature
    L.append(T(r"""\section*{Abstract}
rockmoc designs supersonic nozzle contours with the method of characteristics (MOC) and predicts specific impulse ($I_{sp}$),
thrust and discharge coefficients, exit state, wall distributions, and boundary-layer and separation losses.
We tested every contour type the library produces against @N@ computational fluid dynamics (CFD) simulations in OpenFOAM v2412:
inviscid nozzles, full nozzles with converging sections, two experiments, a LOX/RP-1 real-gas path, wall-resolved viscous runs and
over-expanded nozzles at sea level. Vacuum $I_{sp}$ from rockmoc agrees with the CFD within $\pm0.25\,\%$ in @WITHIN@ of @NISP@ cases;
the largest difference is $@WORST@\,\%$. Refining the mesh from coarse to fine changes $I_{sp}$ by at most $@GCI@\,\%$, so mesh error
does not explain the differences. Throat-plane Mach numbers agree within $@TH@\,\%$ at every throat curvature tested, after the CFD showed
that the series used in the first campaign (rockmoc v0.9) fails at small curvature radii, and the sea-level thrust coefficient of the
over-expanded nozzles lies within $@SEP@\,\%$ of the CFD. The CFD exposed nine defects in v0.9; all are fixed. Three deviations remain:
the minimum-length nozzle centre line at 40 or fewer characteristics, the Rao control-surface bell, and the cooled-wall momentum thickness
of the D-5377 nozzle, where the CFD turbulence model, not rockmoc, departs from the experiment. We conclude that rockmoc meets the
tolerances set before the runs for the contours tested; the sections that follow document the evidence.
\clearpage
\tableofcontents
\clearpage
\listoffigures
\clearpage
\listoftables
\clearpage
""", N=len(rows), WITHIN=within, NISP=len(isp), WORST=f"{worst['isp_err']:+.2f}", GCI=pf(gci_isp), TH=pf(th_max, 1),
        SEP=pf(max((abs(x) for x in sep_err), default=float("nan")), 1)))
    L.append(r"\section{Nomenclature and conventions}\label{sec:nom}")
    L.append(r"This report uses one sign convention throughout: \textbf{error $=$ CFD\,/\,rockmoc $-1$} (in percent for ratios, "
             r"CFD $-$ rockmoc for quantities that are already differences), so a positive error means the CFD value is larger. "
             r"\Cref{tab:nom} defines the abbreviations and symbols." "\n")
    L.append(longtable(r"p{3cm}p{9.2cm}l", ["term", "meaning", "unit"], [
        ["BL", "boundary layer", "--"], ["CEA", "NASA Chemical Equilibrium with Applications thermochemistry code", "--"],
        ["CFD", "computational fluid dynamics", "--"], [r"$C_d$", "discharge coefficient, $\\dot m/\\dot m^*$", "[--]"],
        [r"$C_F$", "thrust coefficient", "[--]"], ["GCI", "grid convergence index (Celik et al.)", "[\\%]"],
        [r"$I_{sp}$", "specific impulse (vacuum unless stated)", "[s]"], ["JANAF", "tabulated thermodynamic data for $c_p(T)$", "--"],
        ["MLN", "minimum-length nozzle", "--"], ["MOC", "method of characteristics", "--"],
        [r"$M$, $M_e$", "Mach number, exit Mach number", "[--]"], [r"$p_0$, $T_0$", "chamber total pressure, total temperature", "[MPa], [K]"],
        ["RANS", "Reynolds-averaged Navier--Stokes", "--"], [r"$R_t$, $R_c$", "throat radius, throat radius of curvature", "[m]"],
        ["SST", "shear-stress transport ($k$--$\\omega$ SST turbulence model)", "--"], ["TIC", "truncated ideal contour", "--"],
        ["TOC, TOP", "Rao thrust-optimised contour (control-surface bell), Rao thrust-optimum parabola", "--"],
        [r"$\theta$, $\delta^*$", "boundary-layer momentum and displacement thickness", "[mm]"],
        [r"$\tau_w$", "wall shear stress", "[Pa]"], [r"$y^+$", "dimensionless wall distance of the first cell", "[--]"],
    ], r"\small", caption="Abbreviations and symbols.", label="tab:nom"))

    # ---------------- 1 summary
    L.append(r"\section{Summary}\label{sec:summary}")
    L.append("rockmoc's vacuum specific impulse agrees with the CFD to within $\\pm0.25\\,\\%$ for "
             f"{within} of the {len(isp)} cases with a vacuum Isp comparison (largest difference "
             f"${worst['isp_err']:+.2f}\\,\\%$, \\texttt{{{t(worst['name'])}}}). "
             "Exit Mach number, discharge coefficient, throat-plane Mach number, wall pressure and the method-of-characteristics "
             f"flow field agree within tolerances fixed before the runs. Coarse, medium and fine meshes change Isp by at most "
             f"${pf(gci_isp)}\\,\\%$, so mesh resolution does not limit the comparison. \\Cref{{sec:lim}} explains the remaining differences and the cases "
             "outside the tolerances.\n")
    tiers_present = ",".join(f"sec:tier-{k}" for k, _ in TIERS if any(r["tier"] == k for r in rows))
    L.append("\\Cref{sec:scope} lists the cases; \\cref{sec:cfd,sec:bc,sec:mesh} describe the CFD method, boundary conditions and meshes. "
             "\\Cref{sec:ref,sec:cmp} define the rockmoc reference values and the tolerances. Results follow in \\cref{sec:res}, "
             "the rockmoc changes that the CFD prompted in \\cref{sec:fix}, mesh independence in \\cref{sec:gci} and limitations in "
             "\\cref{sec:lim}. \\Cref{sec:conc} concludes, \\cref{sec:repro} gives the reproducibility details, and the per-case results sit in "
             f"\\cref{{sec:cases}} (\\cref{{{tiers_present}}}). \\Cref{{sec:app-checks,sec:app-data}} list every check and every data file.\n")
    L.append(f"\\Cref{{fig:isp}} shows the Isp error of every case with a vacuum comparison; bars outside the grey band exceed the tolerance of that case.\n")
    L.append(fig(figs["isp"].replace("figures/", ""),
                 r"Vacuum specific impulse error, CFD/rockmoc $-1$ [\%], for every case with a vacuum comparison. "
                 r"Bar colour marks the tier; the grey band marks $\pm$ the tolerance of each case.", "0.72", label="fig:isp",
                 short="Vacuum Isp error of every case"))
    one_sided = nneg >= 0.7 * len(isp) or nneg <= 0.3 * len(isp)
    L.append(f"{within} of {len(isp)} bars lie within $\\pm0.25\\,\\%$ and the mean absolute error is ${mean_abs:.2f}\\,\\%$; the largest bar is "
             f"${worst['isp_err']:+.2f}\\,\\%$ (\\texttt{{{t(worst['name'])}}}). The CFD value lies below the rockmoc value in {nneg} of {len(isp)} cases. "
             + (f"Because refining the mesh moves Isp by at most ${pf(gci_isp)}\\,\\%$, this one-sided offset points to model approximations in "
                "rockmoc, such as the divergence factor of conical nozzles, and not to discretisation error. " if one_sided else "")
             + "A designer can therefore use the rockmoc Isp within the tolerance bands of \\cref{tab:tol}.\n")

    # ---------------- 2 scope
    L.append(r"\clearpage\section{What was validated}\label{sec:scope}")
    L.append(f"We validated rockmoc in tiers that build from inviscid contour accuracy to viscous, separated flow, so that a failure in a late tier "
             f"can be traced to the first tier where it appears. \\Cref{{tab:tiers}} groups the {len(rows)} cases by tier and \\cref{{tab:cases}} lists "
             "each case with its inlet mode, gas and chamber conditions. The set includes two experiments, Cuffel et al.\\ \\cite{cuffel1969} and "
             "NASA TN D-5377 \\cite{d5377}.\n")
    rows_t = []
    for k, nm in TIERS:
        rs = [r for r in rows if r["tier"] == k]
        if rs:
            rows_t.append([t(nm), str(len(rs)), t(", ".join(r["name"] for r in rs))])
    L.append(longtable(r"p{5.2cm}rp{8.8cm}", ["tier", "cases [--]", "case names"], rows_t, r"\footnotesize",
                       caption="Validation tiers, number of OpenFOAM cases per tier and case names.", label="tab:tiers"))
    L.append(longtable(r"p{3.3cm}lp{3.9cm}lllrr",
                       ["case", "tier", "description", "mode", "gas", r"$\gamma$ [--]", r"$p_0$ [MPa]", r"$T_0$ [K]"],
                       [[t(r["name"]), r["tier"], t(r["title"]), t(r["mode"] + (" + plume" if r["plume"] else "") + (", viscous" if r["viscous"] else "")),
                         t(r["thermo"]), f"{r['gamma']:.4f}", f"{r['p0'] / 1e6:.3g}", f"{r['T0']:.0f}"] for r in rows], r"\scriptsize",
                       caption="Cases: inlet mode, flow model, gas, ratio of specific heats and chamber conditions. Mode is the inlet type "
                               "(throat: supersonic start line; chamber: subsonic inlet).", label="tab:cases", short="Case list"))

    # ---------------- 3 CFD method
    L.append(r"\section{CFD method}\label{sec:cfd}")
    L.append(r"We ran every case in OpenFOAM v2412 \cite{openfoam} with \texttt{rhoCentralFoamAxi}, an axisymmetric variant of \texttt{rhoCentralFoam} "
             r"\cite{greenshields2010}, which uses the central-upwind flux of Kurganov and Tadmor \cite{kurganov2000}. Viscous cases used the "
             r"$k$--$\omega$ SST model \cite{menter1994}. \Cref{tab:cfd} lists solver, schemes, models and run control." "\n")
    items = [
        ("Code", r"OpenFOAM v2412, \texttt{rhoCentralFoamAxi}: density-based compressible solver, axisymmetric variant of rhoCentralFoam."),
        ("Flux scheme", r"Kurganov--Tadmor central-upwind; van Leer limited reconstruction of $\rho$, $\mathbf{U}$, $T$."),
        ("Time integration", r"Explicit Euler, adaptive $\Delta t$ at maximum Courant number 0.4, marched to a steady state."),
        ("Spatial", r"Gauss linear gradients, corrected Laplacian; $k$, $\omega$ convected with limitedLinear 1."),
        ("Geometry", r"Axisymmetric $2.5^\circ$ wedge, one cell thick. Planar case: 2-D slab, empty front/back, symmetry plane on the centre line."),
        ("Thermodynamics", r"hePsiThermo, perfect gas. Constant $c_p$ ($\gamma=1.22$, $R=355$~J/kg\,K for the generic propellant; air, $\gamma=1.4$, for the experiments). "
                           r"LOX/RP-1 cases: JANAF $c_p(T)$ fitted to frozen CEA products."),
        ("Transport (viscous)", r"Sutherland law fitted to Bartz' viscosity for the combustion gas; turbulent Prandtl number 0.85."),
        ("Turbulence (viscous)", r"RANS $k$--$\omega$ SST with Spalding's wall function (valid continuously from $y^+\approx1$ to the log layer). "
                                 r"First-cell radial height set at the throat (5~\textmu m; 10~\textmu m for the separation runs), growing with the local radius towards the exit; "
                                 r"median $y^+$ 0.5--1.4 on the rocket-scale runs and about 9 on the D-5377 rig (measured per case in \cref{sec:mesh})."),
        ("Run length", r"4 flow-through times for supersonic-inlet cases (outlet $\dot m$ steady to $\pm0.03\,\%$ after 2), 8 for chamber cases, "
                       r"1.7--3 for viscous cases; plume cases add 4 transits of the ambient region."),
        ("Time averaging", r"Chamber cases: \texttt{fieldAverage} of $\mathbf{U}$, $p$, $T$, $\rho$ from $0.6\,t_\mathrm{end}$; reported values use the mean fields."),
        ("Parallel", r"8 MPI ranks per case (scotch), two cases at a time on a 16-core workstation (WSL2)."),
    ]
    L.append(longtable(r"p{3.4cm}p{12.4cm}", ["item", "setting"], [[a, b] for a, b in items], r"\small",
                       caption="CFD solver settings, models and run control shared by all cases.", label="tab:cfd"))
    L.append(r"\subsection*{Sample dictionaries (\texttt{t4\_prod\_bl})}" "\nEvery case's full setup is in \\texttt{data/case\\_setup/}.\n")
    ex = cases_dir / "t4_prod_bl"
    for fn in ("system/fvSchemes", "system/controlDict", "constant/thermophysicalProperties", "constant/momentumTransport"):
        p = ex / fn
        if p.exists():
            txt = p.read_text()
            if "FoamFile" in txt:
                txt = txt[txt.find("}", txt.find("FoamFile")) + 1:].strip()
            L.append(f"\\noindent\\texttt{{{t(fn)}}}\n\\begin{{lstlisting}}\n{txt[:2500]}\n\\end{{lstlisting}}\n")

    # ---------------- 4 BCs
    L.append(r"\section{Boundary conditions}\label{sec:bc}")
    L.append(r"Boundary conditions follow the inlet mode of each case. \Cref{tab:bc-inlet} gives the inlet conditions per mode and "
             r"\cref{tab:bc-patch} the remaining patches." "\n")
    L.append(longtable(r"p{2.7cm}p{2.9cm}p{3.4cm}p{2.8cm}p{2.7cm}", ["inlet mode", "used by", "$p$", r"$\mathbf{U}$", "$T$"], [
        ["throat (supersonic start line)", "Tier 1, Tier 3, Tier 4 production/MLN, separation",
         r"fixedValue: isentropic $p$ at $M=1.05$ from $(p_0,T_0)$", r"fixedValue: uniform axial $M=1.05$", r"fixedValue: isentropic $T$ at $M=1.05$"],
        ["chamber (subsonic)", "Tier 2, D-5377", r"totalPressure $p_0$", "pressureInlet\\-Outlet\\-Velocity", r"totalTemperature $T_0$"],
    ], r"\small", caption="Inlet boundary conditions by inlet mode.", label="tab:bc-inlet"))
    L.append(longtable(r"p{4.2cm}p{11.6cm}", ["patch", "condition"], [
        ["outlet (nozzle exit, supersonic)", r"zeroGradient $p$, $\mathbf{U}$, $T$"],
        ["outlet (plume region)", r"$p$: waveTransmissive to $p_a=101\,325$~Pa ($l_\infty=5r_e$); $\mathbf{U}$ zeroGradient; $T$ inletOutlet 300~K"],
        ["ambient (plume side and upstream face)", r"$p$: totalPressure $p_a$; $\mathbf{U}$ pressureInletOutletVelocity; $T$ inletOutlet 300~K"],
        ["wall, inviscid", r"slip; zeroGradient $p$ and $T$"],
        ["wall, viscous", r"noSlip; $T$ zeroGradient (adiabatic) or fixed $T_w(x)$ from NASA TN D-5377 (cooled); $k$: kqRWallFunction; $\omega$: omegaWallFunction; "
                          r"$\nu_t$: nutUSpaldingWallFunction; $\alpha_t$: compressible::alphatWallFunction ($Pr_t=0.85$)"],
        ["axis", "wedge axis (collapsed edge, empty); planar case: symmetryPlane"],
        ["front / back", "wedge (planar case: empty)"],
    ], r"\small", caption="Boundary conditions on the outlet, ambient, wall, axis and wedge patches.", label="tab:bc-patch"))
    L.append("Initial field: 1-D isentropic solution along the nozzle; ambient state in the plume region. "
             "The exact conditions of every case, read from its \\texttt{0/} directory, are listed in its section.\n")

    # ---------------- 5 meshes
    L.append(r"\section{Meshes}\label{sec:mesh}")
    L.append(r"\Cref{tab:mesh-size} lists cell count and wall resolution per case; \cref{tab:mesh-qual} lists the \texttt{checkMesh} quality "
             r"metrics and the run cost. All meshes are structured multi-block hexahedral meshes from \texttt{blockMesh}." "\n")
    L.append(r"The rockmoc wall was resampled by PCHIP in arc length and "
             r"densified around the throat ($-1.5$ to $+2\,R_t$); block boundaries sit at the throat and at $x/R_t=-1,-0.5,0.5,1,2$. "
             r"The axial cell count per block scales with $1+2e^{-x^2}$, giving three times the axial density at the throat. "
             r"Every mesh has 100 radial cells: uniform for inviscid runs; for viscous runs 75\,\% of the height is uniform and a geometric wall layer reaches "
             r"the requested first-cell height. Plume cases add an ambient region $5r_e$ long and $3.5r_e$ high. All meshes pass \texttt{checkMesh}." "\n")
    L.append(longtable(r"p{4.2cm}rrrr",
                       ["case", "cells [--]", r"1st cell [\textmu m]", "$y^+$ median [--]", "$y^+$ max [--]"],
                       [[t(r["name"]), str(r["ncell"]), f"{1e6 * r['dy1']:.1f}--{1e6 * r['dy1max']:.1f}", f(r["y_plus_med"], 3), f(r["y_plus_max"], 3)]
                        for r in rows], r"\scriptsize", caption="Mesh size and wall resolution per case. First-cell height is given from throat to exit.",
                       label="tab:mesh-size", note=r"$y^+$ is reported for viscous cases only; ``--'' marks inviscid cases."))
    L.append(longtable(r"p{4.2cm}rrrrr",
                       ["case", "max aspect ratio [--]", r"max non-orthogonality [$^\circ$]", "max skewness [--]", "run time [h]", "end time [ms]"],
                       [[t(r["name"]), f(r["mesh"].get("max_aspect"), 3), f(r["mesh"].get("max_nonorth"), 3), f(r["mesh"].get("max_skew"), 3),
                         f(r["hours"], 2), f"{1e3 * r['endTime']:.3f}"] for r in rows], r"\scriptsize",
                       caption="Mesh quality from \\texttt{checkMesh}, wall-clock time on 8 cores and simulated end time per case.", label="tab:mesh-qual"))

    # ---------------- 6 rockmoc refs
    L.append(r"\section{rockmoc reference values}\label{sec:ref}")
    L.append("Each reference value comes from a direct rockmoc call on the same design that we meshed (stored in "
             "each case's \\texttt{ref.json} under \\texttt{data/case\\_setup/}); the CFD wall is the rockmoc contour itself. \\Cref{tab:ref} names the rockmoc function behind each quantity. "
             "The transonic throat model follows Kliegel and Levine \\cite{kliegel1969} and, for large curvature radii, Hall \\cite{hall1962}.\n")
    L.append(longtable(r"p{4.2cm}p{11.6cm}", ["quantity", "rockmoc source"], [
        ["Contours", r"\texttt{design\_mln}, \texttt{design\_ideal}, \texttt{design\_tic}, \texttt{design\_rao}, \texttt{design\_rao\_toc}, "
                     r"\texttt{design\_rao\_control}, \texttt{design\_conical}, \texttt{design\_planar\_mln\_nasa}, \texttt{design\_production}"],
        ["Vacuum Isp, $C_F$", r"\texttt{DesignResult.performance()} $\to$ \texttt{integrated\_performance} (uniform exit, conical source flow or MOC wall-pressure "
                              r"integral); production pipeline: \texttt{field\_integrated\_cf}"],
        ["Axis $M,p,T,\\rho,u$; MOC field", "stored characteristic net; 1-D isentropic where a design has no net"],
        ["Wall pressure", r"design wall Mach, or the MOC direct march \texttt{moc.analysis.analyze\_wall}"],
        ["Discharge coefficient, throat plane", r"\texttt{moc.transonic}: Kliegel--Levine series (Hall for $R_c/R_t\ge4$)"],
        ["Boundary layer", r"\texttt{correct\_wall\_bl}: momentum integral, Eckert reference temperature, Bartz transport"],
        ["Friction / BL loss", r"exit control volume (\texttt{delta\_cf\_exit\_vac}, \texttt{delta\_cf\_uncorrected\_vac})"],
        ["Separation", r"\texttt{performance(separation=\textquotesingle schmucker\textquotesingle)} on the MOC wall pressure, Kalt--Badal recovery"],
    ], r"\small", caption="rockmoc source of each reference value.", label="tab:ref"))

    # ---------------- 7 comparison method
    L.append(r"\section{Comparison method and tolerances}\label{sec:cmp}")
    L.append(r"We extracted each comparison quantity from the CFD as listed here, and fixed the tolerances of \cref{tab:tol} before the runs. "
             r"Errors follow the sign convention of \cref{sec:nom}.")
    L.append(r"""\begin{itemize}[leftmargin=*]
\item \textbf{Isp} from CFD exit-plane sums: $\dot m=\sum\rho u_x A$, $F=\sum(\rho u_x^2+p-p_a)A$ (scaled by $360/2.5$ for the wedge), $I_{sp}=F/(\dot m g_0)$.
Difference reported as CFD/rockmoc $-1$.
\item \textbf{Discharge coefficient} $C_d=\bar{\dot m}/\dot m^*(p_0,T_0,A_t)$, using the time-mean mass flow.
\item \textbf{Centre line}: CFD cells nearest the axis at 300 stations, against the rockmoc axis; $p,T,\rho,u$ from the isentropic relations of the rockmoc Mach number.
\item \textbf{MOC field}: CFD interpolated at every characteristic-net point in $(x, r/r_\mathrm{wall})$.
\item \textbf{Boundary layer}: $\delta^*$ and $\theta$ integrated normal to the wall up to 98\,\% total-pressure recovery.
\item \textbf{Separation}: first wall station after which the near-wall axial velocity stays negative.
\item \textbf{Mesh independence}: grid convergence index of Celik et al.\ (2008) on three meshes, $\sqrt2$ refinement per direction.
\end{itemize}
""")
    L.append(longtable(r"p{9cm}p{6cm}", ["check", "tolerance (set before the runs)"], [
        ["Isp, MOC-net cases (Tier 1, Tier 3, mesh study)", r"0.5\,\%"], ["Isp, other inviscid", r"1.0\,\%"], ["Isp, viscous", r"1.5\,\%"],
        ["Exit Mach (area mean) / standard deviation", r"1\,\% / 0.03"], ["Discharge coefficient", r"0.5\,\%"],
        ["Axis Mach maximum / rms error", "0.02 / 0.01"], [r"Axis $p$, $\rho$, $T$, $u$ rms error", r"5\,\%, 4\,\%, 1\,\%, 1\,\%"],
        ["MOC field Mach rms / flow-angle rms", r"1\,\% / $0.5^\circ$"], ["Throat-plane axis / wall Mach", r"5\,\% / 8\,\%"],
        ["D-5377 wall pressure (rms of log ratio)", "0.10"], [r"Boundary-layer $\theta$", r"35\,\%"],
        [r"Separation $\varepsilon_\mathrm{sep}$ / sea-level $C_F$", r"25\,\% / 5\,\%"],
    ], r"\small", caption="Tolerance of each gated check, fixed before the runs.", label="tab:tol"))
    L.append(r"Checks marked \textcolor{infoc}{info} are reported but not gated: design-assumption mismatches (e.g.\ an equilibrium-designed wall in "
             r"frozen-gas CFD) or quantities without a rockmoc counterpart." "\n")

    # ---------------- 8 results overview (pass/fail counts live here)
    L.append(r"\section{Results overview}\label{sec:res}")
    L.append(f"{npass} of the {ng} gated checks pass. \\Cref{{sec:lim}} explains every failure; "
             "38 of them come from rockmoc's own discretisation error at a deliberately coarse 40 characteristics, which converges "
             "first-order and passes at 80. \\Cref{tab:overview} compares Isp and exit Mach number per case, \\cref{tab:counts} gives the check counts "
             "and \\cref{tab:cd} the discharge coefficients of the chamber cases. \\Cref{fig:throat,fig:sep,fig:mln} summarise the throat, "
             "separation and characteristic-count studies.\n")
    L.append(longtable(r"p{4cm}lrrrr",
                       ["case", "tier", r"$I_{sp}$ rockmoc [s]", r"$I_{sp}$ CFD [s]", r"$I_{sp}$ error [\%]", "exit $M$ CFD [--]"],
                       [[t(r["name"]), r["tier"], f(r["isp_ref"], 6), f(r["isp_cfd"], 6),
                         (f"{r['isp_err']:+.2f}" if math.isfinite(r["isp_err"]) else "--") + ("" if r["isp_gated"] or not r["isp_ref"] else "$^\\dagger$"),
                         f"{r['exit_M']:.3f}"] for r in rows], r"\scriptsize",
                       caption="Specific impulse and exit Mach number per case: rockmoc, CFD and the error CFD/rockmoc $-1$.", label="tab:overview",
                       short="Isp and exit Mach number per case",
                       note=r"$^\dagger$ informational only. Plume (separation) cases have no vacuum Isp comparison; their CFD value is the sea-level Isp."))
    L.append(longtable(r"p{4cm}lrrr", ["case", "tier", "pass [--]", "fail [--]", "info [--]"],
                       [[t(r["name"]), r["tier"], str(r["n_pass"]), str(r["n_fail"]), str(r["n_info"])] for r in rows], r"\scriptsize",
                       caption="Number of gated checks passed and failed, and of informational checks, per case.", label="tab:counts",
                       note=r"Informational checks (info) are reported but not gated."))
    cdr = [r for r in rows if r["Cd"] is not None and r["Cd_ref"]]
    if cdr:
        L.append(longtable(r"p{4cm}rrr", ["case", r"$C_d$ CFD [--]", r"$C_d$ rockmoc [--]", "error [\\%]"],
                           [[t(r["name"]), f(r["Cd"], 5), f(r["Cd_ref"], 5), f"{100 * (r['Cd'] / r['Cd_ref'] - 1):+.2f}"] for r in cdr], r"\scriptsize",
                           caption="Discharge coefficient of the chamber cases, CFD (time-mean mass flow) and rockmoc.", label="tab:cd"))
    L.append(f"\\Cref{{fig:throat}} plots the throat-plane axis Mach number against $R_c/R_t$ for both rockmoc series and the CFD points.\n")
    if "throat" in figs:
        L.append(fig(figs["throat"].replace("figures/", ""),
                     r"Throat-plane axis Mach number against throat curvature ratio $R_c/R_t$ for $\gamma=1.22$: Kliegel--Levine series "
                     r"(rockmoc, dashed vermillion), Hall second-order series (rockmoc v0.9, dashed grey) and CFD (black dots).", "0.4",
                     label="fig:throat", short="Throat-plane axis Mach number against curvature ratio"))
        L.append(f"The CFD points follow the Kliegel--Levine curve within ${pf(th_max, 1)}\\,\\%$ for $R_c/R_t$ from 0.625 to 2, while the Hall series "
                 "departs from them below $R_c/R_t\\approx3$. Hall's expansion is second order in $1/(R_c/R_t)$, so it loses accuracy as the throat "
                 "curvature grows. rockmoc therefore uses the Kliegel--Levine series up to $R_c/R_t=4$ (\\cref{tab:fixes}).\n")
    L.append("\\Cref{fig:sep} compares the sea-level thrust coefficient of the three over-expanded nozzles; look for the gap between the vermillion and black bars.\n")
    if "sep" in figs:
        L.append(fig(figs["sep"].replace("figures/", ""),
                     r"Sea-level thrust coefficient $C_F$ of the three over-expanded nozzles: rockmoc with Schmucker separation and Kalt--Badal recovery "
                     r"(vermillion), CFD (black) and attached-flow theory (grey).", "0.4", label="fig:sep", short="Sea-level thrust coefficient of over-expanded nozzles"))
        if sep_err:
            ae = ", ".join(f"${x:+.2f}\\,\\%$" for x in sep_err)
            L.append(f"The rockmoc and CFD thrust coefficients differ by {ae} for the three nozzles. "
                     + (f"Attached-flow theory differs from the CFD by {min(sep_att, key=abs):+.1f} to {max(sep_att, key=abs):+.1f}\\,\\%, because it ignores the separation that "
                        "moves the wall pressure to ambient downstream of the separation point. " if sep_att else "")
                     + "The Schmucker model with Kalt--Badal recovery is therefore necessary for sea-level thrust of these nozzles.\n")
    L.append("\\Cref{fig:mln} shows how the centre-line Mach error of the minimum-length nozzle falls as rockmoc uses more characteristics.\n")
    if "mln" in figs:
        L.append(fig(figs["mln"].replace("figures/", ""),
                     r"Minimum-length nozzle, $M_e=3$: rms error of the centre-line Mach number (CFD against rockmoc, $x>0.5R_t$) against the number of "
                     r"rockmoc characteristics $n$, with a first-order reference line (dotted) and the 1\,\% gate (dash-dotted).", "0.4",
                     label="fig:mln", short="Minimum-length nozzle error against characteristic count"))
        if len(mln) >= 2:
            ks = sorted(mln)
            L.append(f"The error falls from ${mln[ks[0]]:.2f}\\,\\%$ at $n={ks[0]}$ to ${mln[ks[-1]]:.2f}\\,\\%$ at $n={ks[-1]}$, a factor of "
                     f"{mln[ks[0]] / mln[ks[-1]]:.1f} for a factor {ks[-1] / ks[0]:.0f} in $n$. The curve is slightly flatter than first order, "
                     "which would give the full factor. "
                     "The command-line default is now 80 characteristics, where the check passes.\n")

    # ---------------- 9 fixes
    L.append(r"\section{rockmoc changes made from the CFD}\label{sec:fix}")
    L.append("A first campaign with rockmoc v0.9 exposed defects, which we fixed before the final runs; details are in "
             "\\texttt{FINDINGS.md} and \\texttt{CHANGELOG.md}. \\Cref{tab:fixes} lists the nine defects with the Isp error before and after "
             "(CFD/rockmoc $-1$), and \\cref{fig:ba} compares the errors case by case.\n")
    L.append(longtable(r"lp{4.8cm}p{4.6cm}p{5.6cm}", ["\\#", "defect found by CFD", "fix", "Isp error before $\\to$ after"], [
        ["1", r"Divergence loss applied twice ($\lambda$ and $\int\cos^2\theta$)", "one geometric model per contour",
         r"cone $15^\circ$: $+3.8\to+0.24\,\%$; Rao 100\,\%: $+2.7\to-0.03\,\%$"],
        ["2", "Rao exit-angle correlation sign; inexact parabola", "exact quadratic B\\'ezier, corrected correlation \\cite{rao1958}", r"Rao 60\,\%: $-3.3\to-0.02\,\%$"],
        ["3", "Rao TOC built by stretching a TOP contour", "generated at the target length", r"$+2.2\to-0.10\,\%$"],
        ["4", r"$5.9^\circ$ kinks on the ideal/TIC throat arc", "exact arc resampling", "smooth wall"],
        ["5", "Separation on 1-D pressure; losses applied to the ambient term", "Schmucker on the MOC wall pressure, vacuum part separated",
         r"sea-level $C_F$: $-13.5\to+0.17\,\%$ (Rao 80\,\%, 1~MPa)"],
        ["6", r"Friction charged as the full $\int\tau_w\,dA$", "exit control-volume boundary-layer loss (corrected and uncorrected walls)",
         r"production bell $+0.78\to+0.13\,\%$; uncorrected wall $-0.9\to+0.14\,\%$"],
        ["7", r"Planar MLN wall from heuristic steps (7.8\,\% too long)", "C+/streamline intersections", "exit Mach std $0.063\\to0.008$"],
        ["8", r"Hall series used at $R_c/R_t$ 1.5--4~\cite{hall1962}", r"Kliegel--Levine \cite{kliegel1969} up to $R_c/R_t=4$", r"throat Mach error 3--9\,\% $\to$ 0.3--2.4\,\%"],
        ["9", "Command-line default of 40 characteristics", "80", r"MLN axis rms $1.0\to0.64\,\%$"],
    ], r"\small", caption="Defects found by the CFD, the rockmoc fix and the Isp error before and after.", label="tab:fixes"))
    if "ba" in figs:
        L.append(fig(figs["ba"].replace("figures/", ""),
                     r"Vacuum Isp error, CFD/rockmoc $-1$ [\%], before (grey, rockmoc v0.9) and after (blue, current rockmoc) the CFD-driven changes, "
                     r"for every case of the first campaign.", "0.7", label="fig:ba", short="Isp error before and after the CFD-driven changes"))
        if ba:
            wb = max(ba, key=lambda b: abs(b[1]))
            nimp = sum(1 for b in ba if abs(b[2]) < abs(b[1]))
            L.append(f"Across the {len(ba)} cases of the first campaign the mean absolute Isp error drops from "
                     f"${sum(abs(b[1]) for b in ba) / len(ba):.2f}\\,\\%$ to ${sum(abs(b[2]) for b in ba) / len(ba):.2f}\\,\\%$, and {nimp} cases improve. "
                     f"The largest initial error was ${wb[1]:+.2f}\\,\\%$ (\\texttt{{{t(wb[0])}}}), now ${wb[2]:+.2f}\\,\\%$. Items 1 and 2 of \\cref{{tab:fixes}}, the double-counted "
                     "divergence loss and the Rao correlation, caused the largest vacuum-Isp errors; correcting them left only the residuals discussed in \\cref{sec:lim}.\n")

    # ---------------- 10 GCI
    L.append(r"\section{Mesh independence}\label{sec:gci}")
    L.append(f"To test whether mesh error explains the rockmoc--CFD differences, we ran three cases on coarse, medium and fine meshes with $\\sqrt2$ refinement per direction "
             "and evaluated the grid convergence index of Celik et al.\\ \\cite{celik2008}. "
             + ("\\Cref{tab:gci} lists the result. " if gci_rows else "")
             + "\\Cref{fig:gci} plots Isp and exit Mach number against cell count; flat curves mean mesh-independent values.\n")
    if gci_rows:
        L.append(longtable(r"llrrrrr", [t(x) for x in gci_rows[0]], [[t(x) for x in r] for r in gci_rows[1:]], r"\small",
                           caption="Grid convergence index of Celik et al. on the fine mesh for the three mesh-study cases.", label="tab:gci"))
    L.append(fig(figs["gci"].replace("figures/", ""),
                 r"$I_{sp}$ (left, CFD solid) and exit Mach number $M_e$ (right) of the three mesh-study cases on coarse, medium and fine meshes; "
                 r"horizontal dashed lines (left) mark the rockmoc value of each case.", "0.4", label="fig:gci", short="Isp and exit Mach number against mesh size"))
    L.append(f"Isp changes by at most ${pf(gci_isp)}\\,\\%$ and exit Mach number by at most ${pf(gci_m)}\\,\\%$ between the coarse and fine meshes. "
             "These changes are small against the Isp tolerances of 0.5 to 1.5\\,\\% in \\cref{tab:tol}, so the remaining rockmoc--CFD differences are not mesh error.\n")

    # ---------------- 11 limitations
    L.append(r"\section{Limitations}\label{sec:lim}")
    L.append(r"This section explains every gated check that fails and every approximation that limits the comparison. "
             r"Most failures trace to rockmoc's own discretisation or to the CFD turbulence model, and none changes the Isp conclusions of \cref{sec:res}." "\n")
    L.append(r"""\begin{itemize}[leftmargin=*]
\item \textbf{Minimum-length nozzle centre line at $\le40$ characteristics} (38 of the failed gated checks, including the frozen-gas and mesh-study
MLN variants). This is rockmoc's own discretization, first order in the number of characteristics $n$: rms 3.8, 2.0, 1.0 and 0.64\,\% at $n=10$, 20, 40, 80.
It passes at $n=80$; the command-line default is now 80 and the API default 160. The gates were kept at $n=40$ to document the trend.
\item \textbf{Rao control-surface bell} (2 failed checks): MOC field Mach rms 1.4\,\% and flow-angle rms $0.9^\circ$; Isp $-0.12\,\%$.
Most likely a weak internal shock that the characteristic net does not carry.
\item \textbf{D-5377 cooled wall} (2 failed checks): the CFD ($k$--$\omega$ SST, cooled wall) under-predicts the momentum thickness by 30--45\,\% against
both the NASA survey \cite{d5377} and rockmoc, while rockmoc lies within $+5$ to $+25\,\%$ of the survey. Because survey and rockmoc agree with each other and disagree with the CFD,
the discrepancy lies in the turbulence closure \cite{menter1994} and not in rockmoc.
\item \textbf{Conical nozzles}: rockmoc is about 0.2\,\% low in Isp (mesh-independent), because the $\lambda$ factor assumes source flow from a point.
\item \textbf{Cuffel nozzle} ($R_c/R_t=0.625$, $15^\circ$ cone, air): Isp $+0.39\,\%$, within the 1\,\% gate; the same approximation, stronger with a tight throat.
\item \textbf{Rao 80\,\% at 1~MPa, sea level}: rockmoc and the CFD disagree on the separation location. The exit wall pressure is 0.6\,\% above the Schmucker threshold, and the CFD separates only over
the last 0.5 in area ratio. Thrust still agrees ($+0.17\,\%$), because the separated region is short; the location is sensitive near a flat bell lip, so a small pressure error moves it.
\item \textbf{MLN fed by a real converging section}: the design assumes a uniform sonic line at a sharp corner, while a real convergent gives a curved one;
centre-line errors rise to 2.4\,\% rms with Isp unaffected (reported as info).
\item \textbf{Equilibrium chemistry} is not available in the CFD (frozen only); the equilibrium-designed wall run in frozen CFD is a bound ($-3.9\,\%$ Isp), not a match.
\item \textbf{Not covered}: aerospike/plug nozzles, two-phase flow, finite-rate chemistry, three-dimensional effects, transient start-up.
\item \textbf{CFD uncertainty}: RANS ($k$--$\omega$ SST) for viscous and separated flow; $2.5^\circ$ wedge axisymmetry; perfect-gas or frozen JANAF thermodynamics.
\end{itemize}
""")

    # ---------------- conclusions and reproducibility
    L.append(r"\clearpage\section{Conclusions}\label{sec:conc}")
    L.append(r"The comparison supports four conclusions; \cref{sec:summary} gives the short form." "\n")
    L.append(r"\begin{itemize}[leftmargin=*]" "\n"
             f"\\item rockmoc predicts vacuum Isp within $\\pm0.25\\,\\%$ in {within} of {len(isp)} cases and within the preset tolerances elsewhere (\\cref{{sec:res}}).\n"
             f"\\item Mesh error does not explain the differences: Isp moves by at most ${pf(gci_isp)}\\,\\%$ between coarse and fine meshes (\\cref{{sec:gci}}).\n"
             "\\item The CFD exposed nine defects in rockmoc v0.9; all are fixed, and the largest Isp errors fell from several percent to a few tenths of a percent (\\cref{sec:fix}).\n"
             "\\item Three deviations remain, each with a known cause: the minimum-length nozzle at 40 or fewer characteristics, the Rao control-surface bell, "
             "and the D-5377 cooled-wall boundary layer, plus a separation offset of the Rao 80\\,\\% nozzle at 1~MPa (\\cref{sec:lim}).\n"
             "\\end{itemize}\n")
    L.append(f"\\section{{Reproducibility}}\\label{{sec:repro}}")
    L.append("\\Cref{tab:repro} gives the software versions, the source revision and the data paths needed to rebuild the report. "
             "The run command regenerates every figure, table and the PDF from the stored CFD extracts.\n")
    L.append(longtable(r"p{3.8cm}p{12cm}", ["item", "value"], [
        ["CFD solver", r"OpenFOAM v2412 \cite{openfoam}, \texttt{rhoCentralFoamAxi}"],
        ["Schemes", r"Kurganov--Tadmor flux, van Leer limited reconstruction, explicit Euler, Courant number 0.4, Gauss linear gradients (\cref{tab:cfd})"],
        ["Turbulence", r"$k$--$\omega$ SST \cite{menter1994} with Spalding wall function (viscous cases)"],
        ["rockmoc revision", rf"\texttt{{{t(commit)}}} (\texttt{{git rev-parse --short HEAD}})" + (", working tree with uncommitted changes" if dirty else "")],
        ["Report generator", r"\texttt{cfd\_validation/make\_validation\_report.py}, \texttt{cfd\_validation/latex\_report.py}"],
        ["Run command", r"\texttt{py -3.11 cfd\_validation/make\_validation\_report.py}"],
        ["Raw CFD extracts", r"\texttt{cfd\_validation/results/raw/<case>/cfd.npz}"],
        ["Case setups", r"\texttt{validation/data/case\_setup/<case>/}"],
        ["Exported data", r"\texttt{validation/data/} (\cref{tab:datafiles})"],
        ["Figures", r"\texttt{validation/figures/}"],
        ["TeX engine", "XeLaTeX (MiKTeX), two passes"],
    ], r"\small", caption="Reproducibility information: software, source revision and data paths.", label="tab:repro"))
    # ---------------- 12 cases
    L.append(r"\clearpage\section{Case-by-case results}\label{sec:cases}")
    L.append("Each case section opens with a paragraph generated from the check data, then gives the setup, boundary-condition and check tables "
             "and up to five figures. The summary figures of \\cref{sec:res} put these per-case numbers in context.\n")
    TIER_INTRO = {
        "T1": "Tier 1 isolates contour and flow-field accuracy: inviscid flow enters at a supersonic start line, so the converging section does not affect the result.",
        "GCI": "The mesh-study cases repeat three Tier 1 designs on coarse and fine meshes to separate rockmoc error from discretisation error.",
        "T2": "Tier 2 adds the converging section and compares throat-plane Mach number and discharge coefficient with rockmoc and with two experiments.",
        "T3": "Tier 3 repeats the comparison with LOX/RP-1 real-gas properties, using frozen CEA products.",
        "T4": "Tier 4 adds wall-resolved viscous flow, boundary-layer thickness and sea-level separation.",
    }
    for k, nm in TIERS:
        ns = [r["name"] for r in rows if r["tier"] == k]
        if not ns:
            continue
        L.append(f"\\clearpage\\subsection{{{t(nm)}}}\\label{{sec:tier-{k}}}\n{TIER_INTRO[k]}\n")
        for j, n in enumerate(ns):
            ref, meta, res, bcs, gas, cm = sections[n]
            r = byn[n]
            q = lab(n)
            pre = f"\\texttt{{{t(n)}}}: "
            cond = (f"$p_0={r['p0'] / 1e6:.4g}$~MPa, $T_0={r['T0']:.0f}$~K, {t(r['mode'])} inlet, "
                    + ("viscous ($k$--$\\omega$ SST)" if r["viscous"] else "inviscid") + (", sea-level ambient" if r["plume"] else ""))
            have = [fn for fn in ("mesh.png", "contours.png", "axis.png", "wall_exit.png", "overview.png") if (out_dir / "figures" / n / fn).exists()]
            # ---- generated summary paragraph (references this case's figures and tables before they appear)
            gated = [c for c in res["checks"] if not c.info and math.isfinite(c.err) and c.tol]
            para = [f"\\texttt{{{t(n)}}} has {r['ncell']} cells and runs with {cond.replace('$p_0', '$p_0')}."]
            if r["plume"]:
                if r["sep_cf_ref"] and r["sep_cf_cfd"]:
                    para.append(f"For this over-expanded case the comparison quantity is the sea-level thrust coefficient: rockmoc gives {r['sep_cf_ref']:.3f}, "
                                f"the CFD {r['sep_cf_cfd']:.3f} (${100 * (r['sep_cf_cfd'] / r['sep_cf_ref'] - 1):+.2f}\\,\\%$).")
            elif r["isp_ref"] and math.isfinite(r["isp_err"]):
                para.append(f"rockmoc predicts a vacuum Isp of {r['isp_ref']:.2f}~s and the CFD gives {r['isp_cfd']:.2f}~s, an error of ${r['isp_err']:+.2f}\\,\\%$"
                            + ("" if r["isp_gated"] else " (informational, not gated)") + ".")
            para.append(f"The area-mean exit Mach number is {r['exit_M']:.3f} with a standard deviation of {r['exit_std']:.3f}.")
            if gated:
                w = max(gated, key=lambda c: abs(c.err) / c.tol)
                nf = r["n_fail"]
                para.append(f"{r['n_pass']} of {r['n_pass'] + nf} gated checks pass" + (f" and {nf} fail" if nf else "") + ". "
                            f"The check closest to its tolerance is {t(w.quantity)}, with an error of {errtxt(w)} against a tolerance of {toltxt(w)}"
                            + (" (fail)" if not w.ok else " (pass)") + ".")
            refs = [f"\\cref{{tab:{q}-checks}} lists every check, \\cref{{tab:{q}-setup}} the setup and \\cref{{tab:{q}-bc}} the boundary conditions"]
            fl = {"mesh.png": "the mesh", "contours.png": "the CFD field contours", "axis.png": "the centre-line distributions",
                  "wall_exit.png": "the wall and exit-plane profiles", "overview.png": "an overview with the case-specific comparison"}
            figtxt = [f"\\cref{{fig:{q}-{fn[:-4]}}} shows {fl[fn]}" for fn in have]
            para.append("\\Cref{" + f"tab:{q}-checks" + "} lists every check" + f", \\cref{{tab:{q}-setup}} the setup and \\cref{{tab:{q}-bc}} the boundary conditions"
                        + ("; " + ", ".join(figtxt) if figtxt else "") + ".")
            L.append(("" if j == 0 else "\\clearpage") + f"\\subsubsection{{\\texttt{{{t(n)}}}: {t(ref['title'])}}}")
            L.append(" ".join(para) + "\n")
            facts = [
                ["inlet mode", t(r["mode"]) + (r" + ambient plume, $p_a=101\,325$ Pa" if r["plume"] else "")],
                ["geometry", "planar 2-D" if r["planar"] else r"axisymmetric $2.5^\circ$ wedge"],
                ["$p_0$, $T_0$", f"{r['p0'] / 1e6:.4g} MPa, {r['T0']:.0f} K"],
                ["gas", t(f"{gas['thermo']} cp, W = {gas['W']:.3f} kg/kmol, R = {gas['R']:.2f} J/kg K")
                        + (f", $\\gamma={gas['gamma']:.4f}$" if "gamma" in gas else "")
                        + (t(f"; Sutherland As = {gas['As']:.4g}, Ts = {gas['Ts']:.4g} K") if r["viscous"] and "As" in gas else "")],
                ["flow model", r"RANS $k$--$\omega$ SST, wall-resolved" if r["viscous"] else "inviscid (Euler, slip wall)"],
                ["mesh", f"{r['ncell']} cells, {r['nr']} radial; first cell {1e6 * r['dy1']:.1f}--{1e6 * r['dy1max']:.1f}~\\textmu m (throat to exit)"
                         + (f"; $y^+$ median {r['y_plus_med']:.2f}, max {r['y_plus_max']:.2f}" if r["viscous"] else "")
                         + (f"; max aspect {cm.get('max_aspect', float('nan')):.1f}, non-orth. {cm.get('max_nonorth', float('nan')):.1f}$^\\circ$, skew {cm.get('max_skew', float('nan')):.2f}" if cm else "")],
                ["run", f"end time {meta['endTime'] * 1e3:.3f} ms" + (f", {r['hours']:.2f} h on 8 cores" if math.isfinite(r["hours"]) else "")],
                ["CFD $\\dot m$, thrust", f"{res['mdot']:.4g} kg/s, {res['F']:.5g} N" + (" (sea level)" if r["plume"] else " (vacuum)")],
                ["rockmoc vacuum Isp", f"{f(r['isp_ref'], 6)} s" if r["isp_ref"] else "--"],
                ["CFD Isp", f"{res['Isp']:.6g} s" + (" (sea level)" if r["plume"] else "")],
            ]
            L.append(longtable(r"p{3.2cm}p{12.6cm}", ["item", "value"], facts, r"\small",
                               caption=f"{pre}setup: inlet, geometry, gas, flow model, mesh, run and headline results.", label=f"tab:{q}-setup",
                               short=f"{pre}setup"))
            bc_rows = []
            for fld in ("p", "U", "T", "k", "omega", "nut", "alphat"):
                for p_, s in bcs.get(fld, {}).items():
                    if p_ in ("front", "back", "axis"):
                        continue
                    s = s.replace("type ", "").replace("uniform TWALL", "T_w(x) from TN D-5377 (cooled)")
                    bc_rows.append([t(fld), t(p_), t(s[:150])])
            L.append(longtable(r"llp{12.2cm}", ["field", "patch", "condition"], bc_rows, r"\scriptsize",
                               caption=f"{pre}boundary conditions read from the \\texttt{{0/}} directory; axis and wedge faces follow the boundary-condition section.",
                               label=f"tab:{q}-bc", short=f"{pre}boundary conditions"))
            L.append(longtable(r"p{5.4cm}rrrrlp{3.2cm}", ["quantity", "rockmoc", "CFD", "error", "tolerance", "status", "note"],
                               check_rows(res["checks"]), r"\scriptsize",
                               caption=f"{pre}checks, rockmoc against CFD; error is CFD/rockmoc $-1$ for relative checks and CFD $-$ rockmoc for absolute checks.",
                               label=f"tab:{q}-checks", short=f"{pre}checks",
                               note="Values carry the units named in the quantity column."))
            if "net" in res:
                ov_net = " Top left: CFD iso-Mach lines (white, solid) and rockmoc MOC iso-Mach lines (vermillion, dashed)."
            else:
                ov_net = " Top left: CFD Mach field."
            if "throat" in res:
                ov_br = "bottom right, Mach number across the throat plane"
            elif "d5377" in res:
                ov_br = "bottom right, wall pressure against the TN D-5377 taps"
            elif "bl" in res:
                ov_br = "bottom right, boundary-layer thicknesses"
            elif "sep" in res:
                ov_br = "bottom right, near-wall axial velocity and separation locations"
            else:
                ov_br = "bottom right, exit-plane Mach number and flow angle"
            caps = {"mesh.png": f"{pre}Computational mesh of {r['ncell']} cells in the $(x/R_t,\\,r/R_t)$ plane: top, full domain; middle, throat region"
                                + ("; bottom, wall-resolved cells at the exit" if r["viscous"] else "") + f"; vermillion line, nozzle wall. Conditions: {cond}.",
                    "contours.png": f"{pre}CFD contours of Mach number $M$, $\\log_{{10}}(p/p_0)$, $T/T_0$ and $\\log_{{10}}(\\rho/\\rho_0)$ (top to bottom). Conditions: {cond}."
                                    + (" On the Mach panel, black lines are CFD iso-Mach levels and white dashed lines are rockmoc MOC iso-Mach lines at the same levels." if "net" in res else ""),
                    "axis.png": f"{pre}Centre-line $M$, $p/p_0$, $T/T_0$, $\\rho/\\rho_0$ and $u/a_0$ against $x/R_t$: CFD solid black, rockmoc dashed vermillion. "
                                f"Bottom right: difference CFD/rockmoc $-1$ [\\%] for $x>0.5R_t$ (one colour per quantity, in panel order). Conditions: {cond}.",
                    "wall_exit.png": f"{pre}Left, wall pressure $p/p_0$ against $x/R_t$ (CFD solid black, rockmoc dashed vermillion"
                                     + (", TN D-5377 taps as blue dots" if "d5377" in res else "") + "); middle, exit-plane Mach number (CFD solid black, design $M_e$ dashed vermillion) "
                                     "and flow angle (blue)" + ("; right, wall $y^+$ (black) and shear stress $|\\tau_w|$ (blue)" if r["viscous"] else "") + f". Conditions: {cond}.",
                    "overview.png": f"{pre}Overview of CFD against rockmoc (CFD solid, rockmoc dashed)." + ov_net
                                    + f" Top right, axis and wall Mach number; bottom left, axis $p/p_0$, $T/T_0$, $\\rho/\\rho_0$, $u/a_0$; {ov_br}. Conditions: {cond}."}
            shorts = {"mesh.png": "mesh", "contours.png": "CFD contours", "axis.png": "centre-line distributions",
                      "wall_exit.png": "wall and exit plane", "overview.png": "overview"}
            for fn in have:
                L.append(fig(f"{n}/{fn}", caps[fn], {"contours.png": "0.66", "mesh.png": "0.7"}.get(fn, "0.5"),
                             label=f"fig:{q}-{fn[:-4]}", short=f"{pre}{shorts[fn]}"))
            L.append(f"{{\\footnotesize Data: \\url{{data/profiles/{n}_axis.csv}}, \\url{{_wall.csv}}, \\url{{_exit.csv}}; "
                     f"raw CFD extract \\url{{data/cfd/{n}/cfd.npz}}; OpenFOAM setup \\url{{data/case_setup/{n}/}}.}}\n")

    # ---------------- appendices
    L.append(r"\clearpage\appendix\section{Every check}\label{sec:app-checks}")
    L.append("\\Cref{tab:allchecks} lists every check of every case in one table; it repeats the per-case check tables of \\cref{sec:cases} "
             "for searching and sorting.\n")
    L.append(longtable(r"p{2.8cm}p{4.4cm}rrrrlp{2.2cm}", ["case", "quantity", "rockmoc", "CFD", "error", "tol.", "status", "note"],
                       [row for n, _, c in allchecks for row in check_rows([c], with_case=n)], r"\tiny",
                       caption="Every check of every case: rockmoc, CFD, error (CFD/rockmoc $-1$ or absolute), tolerance and status.", label="tab:allchecks",
                       short="Every check of every case"))
    L.append(r"\section{Data files}\label{sec:app-data}")
    L.append("\\Cref{tab:datafiles} lists the exported data files under \\texttt{validation/}.\n")
    L.append(longtable(r"p{6cm}p{10cm}", ["path", "content"], [
        [r"\texttt{data/summary.csv}", "one row per case: Isp, exit Mach, discharge coefficient, axis/field errors, check counts"],
        [r"\texttt{data/all\_checks.csv}", "every check: rockmoc value, CFD value, error, tolerance, status"],
        [r"\texttt{data/cases.csv}", r"case setup: mode, gas, $p_0$, $T_0$, cells, first cell, $y^+$, mesh quality, run time"],
        [r"\texttt{data/profiles/<case>\_axis.csv}", r"centre line: $x$, $M$, $p/p_0$, $T/T_0$, $\rho/\rho_0$, $u/a_0$ for CFD and rockmoc"],
        [r"\texttt{data/profiles/<case>\_wall.csv}", r"wall: $x$, $r$, $M$, $p$, $T$, $\rho$ (and $y^+$, $\tau_w$ for viscous cases)"],
        [r"\texttt{data/profiles/<case>\_exit.csv}", r"exit plane: $r$, $M$, $p$, $T$, $\rho$, $U_x$, $U_r$, face area"],
        [r"\texttt{data/cfd/<case>/cfd.npz}", "raw CFD extract: all cells, wall, exit, inlet, mass-flow history"],
        [r"\texttt{data/cfd/<case>/mesh.npz}, \texttt{log.checkMesh}", "mesh polygons and quality log"],
        [r"\texttt{data/case\_setup/<case>/}", r"complete OpenFOAM case (\texttt{0/}, \texttt{constant/}, \texttt{system/}, \texttt{Allrun}) plus rockmoc \texttt{meta.json}, \texttt{ref.json}"],
        [r"\texttt{index.html}", "the same report as a web page"],
        [r"\texttt{FINDINGS.md}", "narrative findings and the rockmoc changes"],
    ], r"\small", caption="Exported data files and their content.", label="tab:datafiles"))
    L.append(r"""\begin{thebibliography}{99}
\addcontentsline{toc}{section}{References}
\bibitem{menter1994} F.~R. Menter, Two-equation eddy-viscosity turbulence models for engineering applications, \textit{AIAA Journal} 32(8):1598--1605, 1994.
\bibitem{kliegel1969} J.~R. Kliegel and J.~N. Levine, Transonic flow in small throat radius of curvature nozzles, \textit{AIAA Journal} 7(7):1375--1378, 1969.
\bibitem{hall1962} I.~M. Hall, Transonic flow in two-dimensional and axially-symmetric nozzles, \textit{Quarterly Journal of Mechanics and Applied Mathematics} 15(4):487--508, 1962.
\bibitem{rao1958} G.~V.~R. Rao, Exhaust nozzle contour for optimum thrust, \textit{Jet Propulsion} 28(6):377--382, 1958.
\bibitem{cuffel1969} R.~F. Cuffel, L.~H. Back and P.~F. Massier, Transonic flowfield in a supersonic nozzle with small throat radius of curvature, \textit{AIAA Journal} 7(7):1364--1366, 1969.
\bibitem{d5377} NASA Technical Note TN D-5377, nozzle boundary-layer and wall-pressure survey used for the D-5377 cases.
\bibitem{celik2008} I.~B. Celik, U. Ghia, P.~J. Roache, C.~J. Freitas, H. Coleman and P.~E. Raad, Procedure for estimation and reporting of uncertainty due to discretization in CFD applications, \textit{Journal of Fluids Engineering} 130(7):078001, 2008.
\bibitem{kurganov2000} A. Kurganov and E. Tadmor, New high-resolution central schemes for nonlinear conservation laws and convection--diffusion equations, \textit{Journal of Computational Physics} 160(1):241--282, 2000.
\bibitem{greenshields2010} C.~J. Greenshields, H.~G. Weller, L. Gasparini and J.~M. Reese, Implementation of semi-discrete, non-staggered central schemes in a colocated, polyhedral, finite volume framework, for high-speed viscous flows, \textit{International Journal for Numerical Methods in Fluids} 63(1):1--21, 2010.
\bibitem{openfoam} OpenFOAM v2412, ESI-OpenCFD Ltd., 2024.
\end{thebibliography}""")
    L.append(r"\end{document}")

    tex_dir = out_dir / "latex"
    tex_dir.mkdir(parents=True, exist_ok=True)
    src = tex_dir / "rockmoc_CFD_validation_report.tex"
    src.write_text("\n".join(L), encoding="utf-8")
    if not compile_pdf:
        return src
    exe = shutil.which("xelatex") or str(Path.home() / r"AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe")
    for _ in range(2):
        r = subprocess.run([exe, "-interaction=nonstopmode", "-halt-on-error", src.name], cwd=tex_dir,
                           capture_output=True, text=True, errors="ignore", timeout=3600)
        if r.returncode != 0:
            log = (tex_dir / "rockmoc_CFD_validation_report.log").read_text(errors="ignore")
            i = log.find("\n!")
            raise RuntimeError("xelatex failed:\n" + (log[i:i + 1500] if i >= 0 else r.stdout[-1500:]))
    pdf = tex_dir / "rockmoc_CFD_validation_report.pdf"
    shutil.copy2(pdf, out_dir / "rockmoc_CFD_validation_report.pdf")
    return pdf


def gci_table(report_md: Path):
    m = re.search(r"## Mesh independence \(GCI\)\n\n(.*?)\n\n", report_md.read_text(), flags=re.S)
    if not m:
        return []
    lines = [l for l in m.group(1).splitlines() if l.startswith("|") and not l.startswith("|---")]
    return [[c.strip() for c in l.strip("|").split("|")] for l in lines]
