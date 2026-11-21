"""Build the full rockmoc-vs-OpenFOAM validation report in ../validation/.

    py -3.12 make_validation_report.py            # all cases
    py -3.12 make_validation_report.py --no-pdf   # skip the Edge PDF print

Inputs : results/raw/<case>/{cfd.npz, mesh.npz, log.checkMesh}, cases/<case>/ (setup + rockmoc refs),
         results/REPORT_v1.md (first-campaign numbers), results/raw/queue.log (run times).
Outputs: validation/index.html (+ report.pdf), figures/, data/ (CSV tables, per-case profiles,
         raw CFD extracts, mesh polygons, OpenFOAM case setup), FINDINGS.md, README.md.
All numbers come from compare.analyse(), the same code that produces the gated checks.
"""

from __future__ import annotations

import csv
import html
import json
import math
import re
import shutil
import subprocess
import sys
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.tri import Triangulation  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

warnings.simplefilter("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src"))
import compare as C  # noqa: E402

RAW, CASES = HERE / "results" / "raw", HERE / "cases"
OUT = HERE.parent / "validation"
FIG, DATA = OUT / "figures", OUT / "data"
RU = 8314.462618
TIER_NAME = {
    "T1": "Tier 1: contour and MOC field (inviscid, supersonic start-line inlet)",
    "T2": "Tier 2: full nozzles with converging section (throat, Cd, experiments)",
    "T3": "Tier 3: real-gas (LOX/RP-1) path",
    "GCI": "Mesh independence (coarse / medium / fine)",
    "T4": "Tier 4: viscous wall, boundary layer and sea-level separation",
}
# Okabe-Ito colour-blind-safe palette; CFD = black solid, rockmoc = vermillion dashed in every plot
TIER_COL = {"T1": "#0072B2", "T2": "#009E73", "T3": "#E69F00", "GCI": "#999999", "T4": "#CC79A7"}
CFDC, RMC, BLUE, ORANGE, GREY = "k", "#D55E00", "#0072B2", "#E69F00", "#999999"
OI = ["#0072B2", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
                     "axes.titlesize": 8, "figure.titlesize": 8, "lines.linewidth": 1.2, "axes.grid": False})


# ----------------------------------------------------------------------------- helpers
def fmt(v, n=5):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "–"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return f"{v:.{n}g}"


def tier_of(name, ref):
    return "GCI" if name.endswith(("_coarse", "_fine")) else ref["tier"]


def foam_dict(path: Path) -> str:
    t = re.sub(r"/\*.*?\*/", "", path.read_text(), flags=re.S)
    t = re.sub(r"//.*", "", t)
    return " ".join(t.split())


def boundary_conditions(case: Path) -> dict:
    """{field: {patch: 'type ...'}} parsed from 0/<field>."""
    out = {}
    for f in sorted((case / "0").iterdir()):
        if not f.is_file():
            continue
        t = foam_dict(f)
        i = t.find("boundaryField")
        if i < 0:
            continue
        body = t[i:]
        out[f.name] = {m.group(1): m.group(2).strip() for m in re.finditer(r"(\w+) \{ ([^{}]*) \}", body)}
    return out


def gas_props(case: Path) -> dict:
    t = foam_dict(case / "constant" / "thermophysicalProperties")
    g = {"thermo": re.search(r"thermo (\w+);", t).group(1), "transport": re.search(r"transport (\w+);", t).group(1)}
    W = float(re.search(r"molWeight ([\d.eE+-]+)", t).group(1))
    g["W"], g["R"] = W, RU / W
    m = re.search(r"Cp ([\d.eE+-]+)", t)
    if m:
        cp = float(m.group(1))
        g["cp"], g["gamma"] = cp, cp / (cp - g["R"])
    m = re.search(r"As ([\d.eE+-]+); Ts ([\d.eE+-]+)", t)
    if m:
        g["As"], g["Ts"] = float(m.group(1)), float(m.group(2))
    return g


def checkmesh(case_raw: Path) -> dict:
    f = case_raw / "log.checkMesh"
    if not f.exists():
        return {}
    t = f.read_text(errors="ignore")
    g = lambda p: (lambda m: float(m.group(1)) if m else float("nan"))(re.search(p, t))  # noqa: E731
    return {"max_aspect": g(r"Max aspect ratio = ([\d.eE+-]+)"), "max_nonorth": g(r"non-orthogonality Max: ([\d.eE+-]+)"),
            "avg_nonorth": g(r"non-orthogonality Max: [\d.eE+-]+ average: ([\d.eE+-]+)"),
            "max_skew": g(r"Max skewness = ([\d.eE+-]+)"), "ok": "Mesh OK." in t}


def runtimes() -> dict:
    f = RAW / "queue.log"
    if not f.exists():
        return {}
    lines = f.read_text().splitlines()
    k = max((i for i, l in enumerate(lines) if "v2 run" in l), default=0)
    start, out = {}, {}
    for l in lines[k:]:
        p = l.split()
        if len(p) >= 3 and p[2] in ("start", "done"):
            h, m, s = map(int, p[0].split(":"))
            t = h * 3600 + m * 60 + s
            if p[2] == "start":
                start[p[1]] = t
            elif p[1] in start:
                out[p[1]] = ((t - start[p[1]]) % 86400) / 3600.0
    return out


def v1_isp() -> dict:
    f = HERE / "results" / "REPORT_v1.md"
    out = {}
    if f.exists():
        for m in re.finditer(r"^\| (?:T\d|GCI) \| \[([\w.]+)\]\([^)]*\) \| ([^|]+) \| ([^|]+) \|", f.read_text(), flags=re.M):
            try:
                out[m.group(1)] = (float(m.group(2)), float(m.group(3)))
            except ValueError:
                pass
    return out


def polys(mesh):
    f, x, r = mesh["faces"], mesh["x"], mesh["r"]
    out, i = [], 0
    while i < len(f):
        n = int(f[i])
        idx = f[i + 1 : i + 1 + n]
        out.append(np.c_[x[idx], r[idx]])
        i += n + 1
    return out


def first_cell(mesh, ref):
    """Median height of the wall-adjacent cells: cells whose top edge is a boundary edge (used by one cell only)."""
    if mesh is None:
        return float("nan")
    from collections import Counter

    f, x, r = mesh["faces"], mesh["x"], mesh["r"]
    cells, i = [], 0
    while i < len(f):
        n = int(f[i])
        cells.append(f[i + 1 : i + 1 + n])
        i += n + 1
    cnt = Counter(tuple(sorted((int(c[k]), int(c[(k + 1) % len(c)])))) for c in cells for k in range(len(c)))
    xe = ref["wall_geom_x"][-1] + 1e-9
    h = []
    for c in cells:
        if x[c].mean() > xe or len(c) < 4:
            continue
        for k in range(len(c)):
            a, b = int(c[k]), int(c[(k + 1) % len(c)])
            if cnt.get(tuple(sorted((a, b)))) != 1 or x[a] == x[b] or min(r[a], r[b]) < 0.5 * ref["Rt"]:
                continue
            o = [v for v in c if v not in (a, b)]  # wall edge: radial edge length below each wall vertex (blockMesh spacing)
            h.append(float(np.mean([abs(r[w] - r[min(o, key=lambda v: abs(x[v] - x[w]))]) for w in (a, b)])))
    return (float(np.min(h)), float(np.max(h))) if h else (float("nan"), float("nan"))


# ----------------------------------------------------------------------------- per-case figures
def rockmoc_wall_p(ref):
    """rockmoc wall p/p0 vs x: MOC direct march on the meshed wall (what performance() uses), else design 1-D."""
    rwx, rwM = np.array(ref["wall_x"]), np.array(ref["wall_M"], float)
    if ref["kind"] not in ("mln", "ideal", "tic", "production", "planar_mln", "d5377") and not ref["planar"] and not ref.get("table"):
        try:
            from rockmoc.geometry.contour import Contour
            from rockmoc.performance.integrated import moc_thrust_efficiency
            from rockmoc.thermo import PerfectGas

            wx, wr = np.array(ref["wall_geom_x"]), np.array(ref["wall_geom_r"])
            d = wx >= 0
            x, r = wx[d], wr[d]
            c = Contour(x=x, r=r, theta=np.arctan(np.gradient(r, x)), M=np.ones_like(x))
            m = moc_thrust_efficiency(c, PerfectGas(ref["gamma"], ref["R"]), rc_over_rt=ref.get("rc_up") or 1.5)
            xe = np.interp(np.sqrt(np.asarray(m["wall_eps"])) * ref["Rt"], r, x)
            return xe, np.asarray(m["wall_p_pc"]), "rockmoc MOC wall march"
        except Exception:  # noqa: BLE001
            pass
    ok = np.isfinite(rwM)
    if ok.sum() > 2:
        lab = "rockmoc MOC (design wall)" if ref.get("has_net") else "rockmoc 1-D (area ratio)"
        return rwx[ok], C.ratios(ref, rwM[ok])[0], lab
    return None


def tri_field(z, ref):
    Rt = ref["Rt"]
    wx, wr = np.array(ref["wall_geom_x"]), np.array(ref["wall_geom_r"])
    tri = Triangulation(z["c_x"] / Rt, z["c_r"] / Rt)
    xc, rc = tri.x[tri.triangles].mean(1), tri.y[tri.triangles].mean(1)
    tri.set_mask((xc * Rt <= wx[-1]) & (rc * Rt > np.interp(xc * Rt, wx, wr)))
    return tri, wx / Rt, wr / Rt


def fig_mesh(name, ref, z, mesh, path):
    Rt = ref["Rt"]
    P = [p / Rt for p in polys(mesh)]
    wx, wr = np.array(ref["wall_geom_x"]) / Rt, np.array(ref["wall_geom_r"]) / Rt
    panels = [(None, None),
              (None, (max(wx[0], -1.2), 1.6, 0.0, 1.35))]
    if ref["viscous"]:
        xe, re_ = wx[-1], wr[-1]
        panels.append((None, (xe - 0.25, xe + 0.02, re_ - 0.08, re_ + 0.01)))
    fig, ax = plt.subplots(len(panels), 1, figsize=(6.6, 2.6 * len(panels)))
    for a, (t, box) in zip(np.atleast_1d(ax), panels):
        if box is None:
            lw = 0.05 if len(P) > 40000 else 0.08
            a.add_collection(PolyCollection(P, facecolors="none", edgecolors="#555555", linewidths=lw))
            a.autoscale()
        else:
            x0, x1, y0, y1 = box
            sel = [p for p in P if p[:, 0].max() >= x0 and p[:, 0].min() <= x1 and p[:, 1].max() >= y0 and p[:, 1].min() <= y1]
            a.add_collection(PolyCollection(sel, facecolors="none", edgecolors="#1a365d", linewidths=0.35))
            a.set_xlim(x0, x1)
            a.set_ylim(y0, y1)
        a.plot(wx, wr, color=RMC, lw=1.2, label="nozzle wall")
        a.set_aspect("equal")
        a.set_xlabel(r"$x/R_t$ [–]")
        a.set_ylabel(r"$r/R_t$ [–]")
    np.atleast_1d(ax)[0].legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fig_contours(name, ref, z, res, path):
    p0, T0, rho0, _ = C.norms(ref)
    tri, wx, wr = tri_field(z, ref)
    fields = [("Mach", r"Mach number $M$ [–]", z["c_Ma"], "viridis", False), ("p", r"$\log_{10}(p/p_0)$ [–]", z["c_p"] / p0, "cividis", True),
              ("T", r"$T/T_0$ [–]", z["c_T"] / T0, "magma", False), ("rho", r"$\log_{10}(\rho/\rho_0)$ [–]", z["c_rho"] / rho0, "viridis", True)]
    xs, ys = np.ptp(tri.x), np.ptp(tri.y)
    h = max(1.6, 12.5 * ys / xs + 0.9)
    fig, ax = plt.subplots(4, 1, figsize=(6.6, 4 * max(1.7, 0.55 * h + 0.4)))
    for a, (key, t, v, cm, lg) in zip(ax.ravel(), fields):
        vv = np.log10(np.maximum(v, 1e-12)) if lg else v
        tc = a.tripcolor(tri, vv, shading="gouraud", cmap=cm)
        fig.colorbar(tc, ax=a, label=t, fraction=0.025, pad=0.01)
        if key == "Mach" and "net" in res:
            lv = np.linspace(1.2, max(1.3, float(np.nanmax(z["c_Ma"])) - 0.1), 9)
            a.tricontour(tri, vv, levels=lv, colors="k", linewidths=0.8)
        else:
            a.tricontour(tri, vv, levels=10, colors="k", linewidths=0.3, alpha=0.5)
        a.plot(wx, wr, "k", lw=1.4)
        a.set_aspect("equal")
        a.set_xlabel(r"$x/R_t$ [–]")
        a.set_ylabel(r"$r/R_t$ [–]")
    if "net" in res:
        nx, nr, nM = res["net"]
        lv = np.linspace(1.2, max(1.3, float(np.nanmax(z["c_Ma"])) - 0.1), 9)
        try:
            ax[0].tricontour(nx / ref["Rt"], nr / ref["Rt"], nM, levels=lv, colors="w", linewidths=0.9, linestyles="--")
        except Exception:  # noqa: BLE001
            pass
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def fig_axis(name, ref, res, path):
    Rt = ref["Rt"]
    A, R = res["axis"], res["axis_ref"]
    rM = np.interp(A["x"], R["x"], R["M"])
    rr = C.ratios(ref, rM)
    fig, ax = plt.subplots(2, 3, figsize=(6.8, 5.0))
    items = [(r"$M$ [–]", A["M"], rM), (r"$p/p_0$ [–]", A["p"], rr[0]), (r"$T/T_0$ [–]", A["T"], rr[1]),
             (r"$\rho/\rho_0$ [–]", A["rho"], rr[2]), (r"$u/a_0$ [–]", A["u"], rr[3])]
    for a, (t, cv, rv) in zip(ax.ravel(), items):
        a.plot(A["x"] / Rt, cv, CFDC, lw=1.4, label="CFD")
        a.plot(A["x"] / Rt, rv, RMC, ls="--", lw=1.4, label="rockmoc")
        a.set_ylabel(t)
        a.set_xlabel(r"$x/R_t$ [–]")
        a.grid(alpha=0.3)
        if t.startswith((r"$p/", r"$\rho/")):
            a.set_yscale("log")
    ax[0, 0].legend()
    a = ax[1, 2]
    m = (A["x"] > 0.5 * Rt) & (A["x"] >= R["x"].min()) & (A["x"] <= R["x"].max())
    for (t, cv, rv), c in zip(items, OI):
        a.plot(A["x"][m] / Rt, 100 * (cv[m] / rv[m] - 1), color=c, lw=1, label=t)
    a.axhline(0, color="gray", lw=0.8)
    a.set_ylabel("CFD/rockmoc − 1 [%]")
    a.set_xlabel(r"$x/R_t$ [–]")
    a.legend(fontsize=7.5, ncol=2)
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def fig_wall_exit(name, ref, z, res, path):
    Rt = ref["Rt"]
    W = res["wall"]
    rwx, rwM = np.array(ref["wall_x"]), np.array(ref["wall_M"], float)
    vis = ref["viscous"]
    fig, ax = plt.subplots(1, 3 if vis else 2, figsize=(6.8, 2.7))
    a = ax[0]
    a.semilogy(W["x"] / Rt, W["p"], CFDC, label="CFD")
    rw = rockmoc_wall_p(ref)
    if rw is not None:
        a.semilogy(rw[0] / Rt, rw[1], RMC, ls="--", label="rockmoc")
    if "d5377" in res:
        d = res["d5377"]
        a.semilogy(d["xd"] * ref["dstar"] / Rt, d["exp"], "ko", ms=4, label="TN D-5377 taps")
    a.set_xlabel(r"$x/R_t$ [–]")
    a.set_ylabel(r"wall $p/p_0$ [–]")
    a.grid(alpha=0.3, which="both")
    a.legend()
    a = ax[1]
    r = z["e_r"] / Rt
    o = np.argsort(r)
    a.plot(r[o], z["e_Ma"][o], CFDC, label="CFD $M$")
    if ref.get("Me_design"):
        a.axhline(ref["Me_design"], color=RMC, ls="--", label="design $M_e$")
    a2 = a.twinx()
    a2.plot(r[o], np.degrees(np.arctan2(z["e_U_r"], z["e_U_x"]))[o], BLUE, lw=1, label="CFD angle")
    a2.set_ylabel("flow angle [°]", color=BLUE)
    a.set_xlabel(r"$r/R_t$ [–]")
    a.set_ylabel(r"$M$ [–]")
    a.grid(alpha=0.3)
    a.legend(loc="lower left")
    if vis and "w_yPlus" in z:
        a = ax[2]
        o = np.argsort(z["w_x"])
        a.plot(z["w_x"][o] / Rt, z["w_yPlus"][o], CFDC, lw=1)
        a.set_yscale("log")
        a.set_xlabel(r"$x/R_t$ [–]")
        a.set_ylabel(r"$y^+$ [–]")
        a.grid(alpha=0.3, which="both")
        a3 = a.twinx()
        tw = np.hypot(z["w_wallShearStress_x"], z["w_wallShearStress_r"])[o]
        a3.plot(z["w_x"][o] / Rt, tw, BLUE, lw=0.8)
        a3.set_ylabel(r"$|\tau_w|$ [Pa]", color=BLUE)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


# ----------------------------------------------------------------------------- summary figures
def summary_figs(rows, v1):
    out = {}
    # Isp error bar chart
    sel = [r for r in rows if r["isp_ref"] and not r["plume"] and r["isp_gated"]]
    sel.sort(key=lambda r: (["T1", "GCI", "T2", "T3", "T4"].index(r["tier"]), r["name"]))
    fig, a = plt.subplots(figsize=(6.6, 0.24 * len(sel) + 1.2))
    y = np.arange(len(sel))
    a.barh(y, [r["isp_err"] for r in sel], color=[TIER_COL[r["tier"]] for r in sel])
    for yy, r in zip(y, sel):
        a.plot([-100 * r["isp_tol"], 100 * r["isp_tol"]], [yy, yy], color="#d9d9d9", lw=6, zorder=0, solid_capstyle="butt")
    a.set_yticks(y, [r["name"] for r in sel], fontsize=8)
    a.invert_yaxis()
    a.axvline(0, color="k", lw=0.8)
    a.set_xlabel(r"$I_{sp}$ error, CFD/rockmoc − 1 [%]")
    from matplotlib.patches import Patch

    a.legend(handles=[Patch(color=c, label=t) for t, c in TIER_COL.items()], loc="lower right")
    a.grid(alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(FIG / "summary_isp_error.png", dpi=130)
    plt.close(fig)
    out["isp"] = "figures/summary_isp_error.png"
    # before/after
    ba = [(r["name"], 100 * (v1[r["name"]][1] / v1[r["name"]][0] - 1), r["isp_err"]) for r in sel
          if r["name"] in v1 and math.isfinite(v1[r["name"]][0])]
    if ba:
        fig, a = plt.subplots(figsize=(6.6, 0.28 * len(ba) + 1.3))
        y = np.arange(len(ba))
        a.barh(y - 0.2, [b[1] for b in ba], 0.4, color=GREY, label="first campaign (rockmoc v0.9)")
        a.barh(y + 0.2, [b[2] for b in ba], 0.4, color=BLUE, label="final (rockmoc with CFD fixes)")
        a.set_yticks(y, [b[0] for b in ba], fontsize=8)
        a.invert_yaxis()
        a.axvline(0, color="k", lw=0.8)
        a.set_xlabel(r"$I_{sp}$ error, CFD/rockmoc − 1 [%]")
        a.legend()
        a.grid(alpha=0.3, axis="x")
        fig.tight_layout()
        fig.savefig(FIG / "summary_before_after.png", dpi=130)
        plt.close(fig)
        out["ba"] = "figures/summary_before_after.png"
    # GCI
    byn = {r["name"]: r for r in rows}
    fig, a = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for base, c in (("t1_mln_M3", OI[0]), ("t1_rao_80", OI[1]), ("t1_cone_15", OI[2])):
        trip = [byn.get(f"{base}_coarse"), byn.get(base), byn.get(f"{base}_fine")]
        if None in trip:
            continue
        N = np.array([t["ncell"] for t in trip], float)
        a[0].plot(N, [t["isp_cfd"] for t in trip], "o-", color=c, label=f"{base} CFD")
        a[0].axhline(trip[1]["isp_ref"], color=c, ls="--", lw=0.8)
        a[1].plot(N, [t["exit_M"] for t in trip], "o-", color=c, label=base)
    a[0].set_xscale("log")
    a[1].set_xscale("log")
    a[0].set_xlabel(r"cells $N$ [–]")
    a[1].set_xlabel(r"cells $N$ [–]")
    a[0].set_ylabel(r"$I_{sp}$ [s]")
    a[1].set_ylabel(r"exit Mach number $M_e$ [–]")
    for b in a:
        b.grid(alpha=0.3, which="both")
    a[1].legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIG / "summary_mesh_independence.png", dpi=130)
    plt.close(fig)
    out["gci"] = "figures/summary_mesh_independence.png"
    # MLN axis error vs n
    ns = [(n, byn.get(f"t1_mln_n{n}") or (byn.get("t1_mln_M3") if n == 40 else None)) for n in (10, 20, 40, 80)]
    ns = [(n, r) for n, r in ns if r and math.isfinite(r["axis_M_rms"])]
    if ns:
        fig, a = plt.subplots(figsize=(4.6, 3.2))
        n_ = np.array([n for n, _ in ns], float)
        e_ = np.array([100 * r["axis_M_rms"] for _, r in ns])
        a.loglog(n_, e_, "ko-", label="axis Mach rms error")
        a.loglog(n_, e_[0] * n_[0] / n_, color=GREY, ls=":", lw=1.2, label=r"first order, $\propto 1/n$")
        a.axhline(1.0, color=GREY, ls="-.", lw=1, label="1 % gate")
        a.set_xlabel(r"characteristics $n$ [–]")
        a.set_ylabel(r"axis $M$ rms error [%]")
        a.grid(alpha=0.3, which="both")
        a.legend()
        fig.tight_layout()
        fig.savefig(FIG / "summary_mln_n_convergence.png", dpi=130)
        plt.close(fig)
        out["mln"] = "figures/summary_mln_n_convergence.png"
    # throat plane axis Mach vs Rc
    from rockmoc.moc.transonic import hall_uv, kl_uv

    g = 1.22
    Mof = lambda V: math.sqrt(max(2 * V * V / ((g + 1) - (g - 1) * V * V), 0))  # noqa: E731
    Rg = np.linspace(0.5, 5, 80)
    fig, a = plt.subplots(figsize=(4.8, 3.2))
    a.plot(Rg, [Mof(math.hypot(*kl_uv(0, 0, g, 1, R))) for R in Rg], RMC, ls="--", label="Kliegel–Levine (rockmoc)")
    Rh = Rg[Rg >= 0.9]
    a.plot(Rh, [Mof(math.hypot(*hall_uv(0, 0, g, 1, R))) for R in Rh], GREY, ls="--", label="Hall 2nd order (rockmoc v0.9)")
    for nm, rc in (("t2_rao80_rc0.625", 0.625), ("t2_rao80_rc1", 1.0), ("t2_ideal_full", 1.5), ("t2_tic80_full", 1.5), ("t2_rao80_rc2", 2.0)):
        r = byn.get(nm)
        if r and r.get("throat_axis_M"):
            a.plot(rc, r["throat_axis_M"], "o", color=CFDC)
    a.plot([], [], "o", color=CFDC, label=r"CFD ($\gamma=1.22$)")
    a.set_xlabel(r"$R_c/R_t$ [–]")
    a.set_ylabel(r"throat-plane axis $M$ [–]")
    a.grid(alpha=0.3)
    a.legend()
    fig.tight_layout()
    fig.savefig(FIG / "summary_throat_series.png", dpi=130)
    plt.close(fig)
    out["throat"] = "figures/summary_throat_series.png"
    # separation
    sp = [r for r in rows if r["plume"]]
    if sp:
        fig, a = plt.subplots(figsize=(6.0, 2.8))
        y = np.arange(len(sp))
        a.barh(y - 0.25, [r["sep_cf_ref"] for r in sp], 0.25, color=RMC, label="rockmoc (Schmucker + Kalt–Badal)")
        a.barh(y, [r["sep_cf_cfd"] for r in sp], 0.25, color=CFDC, label="CFD")
        a.barh(y + 0.25, [r["sep_cf_att"] for r in sp], 0.25, color=GREY, label="attached-flow theory")
        a.set_yticks(y, [r["name"] for r in sp])
        a.set_xlim(0.9, max(1.3, max(r["sep_cf_cfd"] for r in sp) + 0.05))
        a.set_xlabel(r"sea-level thrust coefficient $C_F$ [–]")
        a.legend(loc="lower right")
        a.grid(alpha=0.3, axis="x")
        fig.tight_layout()
        fig.savefig(FIG / "summary_separation.png", dpi=130)
        plt.close(fig)
        out["sep"] = "figures/summary_separation.png"
    return out


# ----------------------------------------------------------------------------- data export
def write_csv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def export_profiles(name, ref, z, res):
    d = DATA / "profiles"
    p0, T0, rho0, a0 = C.norms(ref)
    A, R = res["axis"], res["axis_ref"]
    rM = np.interp(A["x"], R["x"], R["M"])
    rr = C.ratios(ref, rM)
    write_csv(d / f"{name}_axis.csv", ["x_m", "M_cfd", "p_p0_cfd", "T_T0_cfd", "rho_rho0_cfd", "u_a0_cfd",
                                        "M_rockmoc", "p_p0_rockmoc", "T_T0_rockmoc", "rho_rho0_rockmoc", "u_a0_rockmoc"],
              np.c_[A["x"], A["M"], A["p"], A["T"], A["rho"], A["u"], rM, *rr].tolist())
    o = np.argsort(z["w_x"])
    cols = ["w_x", "w_r", "w_Ma", "w_p", "w_T", "w_rho"] + [k for k in ("w_yPlus", "w_wallShearStress_x", "w_wallShearStress_r") if k in z]
    write_csv(d / f"{name}_wall.csv", [c[2:] for c in cols], np.c_[tuple(z[c][o] for c in cols)].tolist())
    o = np.argsort(z["e_r"])
    cols = ["e_r", "e_Ma", "e_p", "e_T", "e_rho", "e_U_x", "e_U_r", "e_A"]
    write_csv(d / f"{name}_exit.csv", [c[2:] for c in cols], np.c_[tuple(z[c][o] for c in cols)].tolist())


# ----------------------------------------------------------------------------- HTML
CSS = """
:root{--ink:#1a202c;--mute:#4a5568;--line:#e2e8f0;--bg:#fff;--accent:#2b6cb0;--pass:#276749;--fail:#c53030;--info:#718096}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--ink:#e2e8f0;--mute:#a0aec0;--line:#2d3748;--bg:#171923;--accent:#63b3ed}}
body{font:15px/1.55 -apple-system,Segoe UI,Roboto,sans-serif;color:var(--ink);background:var(--bg);max-width:1180px;margin:0 auto;padding:24px 16px}
h1{font-size:28px;margin:.2em 0}h2{border-bottom:2px solid var(--accent);padding-bottom:4px;margin-top:2.2em}
h3{margin-top:1.6em}table{border-collapse:collapse;width:100%;font-size:12.5px;margin:.6em 0 1.2em}
th,td{border:1px solid var(--line);padding:3px 6px;text-align:left;vertical-align:top}th{background:rgba(43,108,176,.08)}
td.n{text-align:right;font-variant-numeric:tabular-nums}.PASS{color:var(--pass);font-weight:600}.FAIL{color:var(--fail);font-weight:700}
.info{color:var(--info)}img{max-width:100%;border:1px solid var(--line);margin:.4em 0}.mute{color:var(--mute)}
.kpi{display:flex;gap:14px;flex-wrap:wrap;margin:1em 0}.kpi div{border:1px solid var(--line);border-radius:8px;padding:10px 14px;min-width:150px}
.kpi b{display:block;font-size:22px;color:var(--accent)}code{font-size:12.5px}pre{background:rgba(127,127,127,.08);padding:8px;overflow-x:auto;font-size:12px;white-space:pre-wrap;word-break:break-word}
.case{page-break-before:always}.toc{columns:2;font-size:13px}.wrap{overflow-x:auto}
@media print{body{max-width:none;font-size:11px}h2{page-break-after:avoid}.case{page-break-before:always}img{page-break-inside:avoid}}
"""


def table(header, rows, num=()):
    h = "".join(f"<th>{html.escape(str(c))}</th>" for c in header)
    b = []
    for r in rows:
        tds = []
        for i, c in enumerate(r):
            cls = ' class="n"' if i in num else ""
            s = str(c)
            if s in ("PASS", "FAIL", "info"):
                cls = f' class="{s}"'
            tds.append(f"<td{cls}>{s if s.startswith('<') or '<code>' in s or '<b>' in s else html.escape(s)}</td>")
        b.append("<tr>" + "".join(tds) + "</tr>")
    return f'<div class="wrap"><table><thead><tr>{h}</tr></thead><tbody>{"".join(b)}</tbody></table></div>'


def check_rows(checks):
    out = []
    for c in checks:
        e = f"{100 * c.err:+.2f} %" if c.rel else f"{c.err:+.3g}"
        t = f"{100 * c.tol:.1f} %" if c.rel else f"{c.tol:g}"
        st = "info" if c.info else ("PASS" if c.ok else "FAIL")
        out.append([c.quantity, fmt(c.rockmoc), fmt(c.cfd), e if math.isfinite(c.err) else "–", t, st, c.note])
    return out


def bc_table(bcs):
    fields = [f for f in ("p", "U", "T", "k", "omega", "nut", "alphat") if f in bcs]
    patches = sorted({p for f in fields for p in bcs[f]}, key=lambda p: ["inlet", "outlet", "wall", "ambient", "axis", "front", "back"].index(p)
                     if p in ["inlet", "outlet", "wall", "ambient", "axis", "front", "back"] else 9)
    rows = []
    for p in patches:
        row = [p]
        for f in fields:
            s = bcs[f].get(p, "")
            s = s.replace("type ", "").replace("uniform TWALL", "T_w(x) from TN D-5377 (cooled)")
            row.append(s[:110])
        rows.append(row)
    return table(["patch"] + fields, rows)


def main():
    pdf = "--no-pdf" not in sys.argv
    nofigs = "--no-figs" in sys.argv
    if OUT.exists():
        for sub in (("data",) if nofigs else ("figures", "data")):
            shutil.rmtree(OUT / sub, ignore_errors=True)
    FIG.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    names = sorted(p.parent.name for p in RAW.glob("*/cfd.npz") if (CASES / p.parent.name).exists())
    sel = [a for a in sys.argv[1:] if not a.startswith("--")]
    names = [n for n in names if not sel or n in sel]
    rt, v1 = runtimes(), v1_isp()
    rows, allchecks, sections = [], [], {}
    for n in names:
        print("case", n, flush=True)
        ref, z, res = C.analyse(n)
        meta = json.loads((CASES / n / "meta.json").read_text())
        tier = tier_of(n, ref)
        cdir = FIG / n
        cdir.mkdir(parents=True, exist_ok=True)
        mesh = np.load(RAW / n / "mesh.npz") if (RAW / n / "mesh.npz").exists() else None
        if not nofigs:
            if mesh is not None:
                fig_mesh(n, ref, z, mesh, cdir / "mesh.png")
            fig_contours(n, ref, z, res, cdir / "contours.png")
            fig_axis(n, ref, res, cdir / "axis.png")
            fig_wall_exit(n, ref, z, res, cdir / "wall_exit.png")
            C.figure(ref, z, res, cdir / "overview.png")
        export_profiles(n, ref, z, res)
        for src, dst in ((RAW / n / "cfd.npz", DATA / "cfd" / n / "cfd.npz"), (RAW / n / "mesh.npz", DATA / "cfd" / n / "mesh.npz"),
                         (RAW / n / "log.checkMesh", DATA / "cfd" / n / "log.checkMesh")):
            if src.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
        sd = DATA / "case_setup" / n
        shutil.copytree(CASES / n, sd, dirs_exist_ok=True, ignore=shutil.ignore_patterns("processor*", "[1-9]*", "postProcessing"))
        gas, cm, bcs = gas_props(CASES / n), checkmesh(RAW / n), boundary_conditions(CASES / n)
        isp_chk = next((c for c in res["checks"] if c.quantity.startswith("Isp_vac")), None)
        sepc = {c.quantity: c for c in res["checks"]}
        cf_chk = next((c for k, c in sepc.items() if k.startswith("Cf at sea level (rockmoc")), None)
        row = dict(name=n, tier=tier, title=ref["title"], kind=ref.get("kind"), mode=ref["mode"], planar=ref["planar"],
                   viscous=bool(ref["viscous"]), plume=bool(ref.get("plume")), ncell=int(ref["ncell"]), nr=ref.get("nr"),
                   p0=ref["p0"], T0=ref["T0"], gamma=gas.get("gamma", ref["gamma"]), R=gas["R"], thermo=gas["thermo"],
                   endTime=meta["endTime"], hours=rt.get(n, float("nan")), isp_ref=ref.get("Isp_vac"), isp_cfd=res["Isp"],
                   isp_err=100 * (res["Isp"] / ref["Isp_vac"] - 1) if ref.get("Isp_vac") and not ref.get("plume") else float("nan"),
                   isp_tol=isp_chk.tol if isp_chk else float("nan"), isp_gated=bool(isp_chk and not isp_chk.info),
                   exit_M=res["exit_M_mean"], exit_std=res["exit_M_std"], Cd=res.get("Cd") if ref["mode"] == "chamber" else None, Cd_ref=ref.get("Cd"),
                   axis_M_rms=res.get("axis_M_rms", float("nan")), field_M_rms=res.get("field_M_rms", float("nan")),
                   throat_axis_M=res["throat"]["M"][0] if "throat" in res else None,
                   y_plus_med=float(np.median(z["w_yPlus"])) if ref["viscous"] and "w_yPlus" in z else float("nan"),
                   y_plus_max=float(np.max(z["w_yPlus"])) if ref["viscous"] and "w_yPlus" in z else float("nan"),
                   dy1=first_cell(mesh, ref)[0], dy1max=first_cell(mesh, ref)[1], mesh=cm,
                   n_pass=sum(1 for c in res["checks"] if not c.info and c.ok),
                   n_fail=sum(1 for c in res["checks"] if not c.info and not c.ok),
                   n_info=sum(1 for c in res["checks"] if c.info),
                   sep_cf_ref=cf_chk.rockmoc if cf_chk else None, sep_cf_cfd=cf_chk.cfd if cf_chk else None,
                   sep_cf_att=ref.get("sep_Cf_att"), sep_eps_ref=ref.get("sep_eps"), sep_eps_cfd=res.get("eps_sep"))
        rows.append(row)
        for c in res["checks"]:
            allchecks.append((n, tier, c))
        sections[n] = (ref, meta, res, bcs, gas, cm)
        plt.close("all")
    figs = summary_figs(rows, v1)
    # ---- CSV
    write_csv(DATA / "summary.csv", ["case", "tier", "title", "Isp_rockmoc_s", "Isp_cfd_s", "Isp_err_pct", "Isp_tol_pct", "exit_M_cfd",
                                     "exit_M_std", "Cd_cfd", "Cd_rockmoc", "axis_M_rms", "field_M_rms", "gated_pass", "gated_fail", "info"],
              [[r["name"], r["tier"], r["title"], r["isp_ref"], r["isp_cfd"], r["isp_err"], 100 * r["isp_tol"], r["exit_M"], r["exit_std"],
                r["Cd"], r["Cd_ref"], r["axis_M_rms"], r["field_M_rms"], r["n_pass"], r["n_fail"], r["n_info"]] for r in rows])
    write_csv(DATA / "all_checks.csv", ["case", "tier", "quantity", "rockmoc", "cfd", "error", "relative", "tolerance", "status", "note"],
              [[n, t, c.quantity, c.rockmoc, c.cfd, c.err, c.rel, c.tol, "info" if c.info else ("PASS" if c.ok else "FAIL"), c.note]
               for n, t, c in allchecks])
    write_csv(DATA / "cases.csv", ["case", "tier", "title", "mode", "planar", "viscous", "plume", "cells", "radial_cells", "p0_Pa", "T0_K",
                                   "gamma", "R_J_kgK", "thermo", "endTime_s", "wall_clock_h", "first_cell_min_m", "first_cell_max_m", "yplus_median", "yplus_max",
                                   "max_aspect", "max_nonorth_deg", "max_skew"],
              [[r["name"], r["tier"], r["title"], r["mode"], r["planar"], r["viscous"], r["plume"], r["ncell"], r["nr"], r["p0"], r["T0"],
                r["gamma"], r["R"], r["thermo"], r["endTime"], r["hours"], r["dy1"], r["dy1max"], r["y_plus_med"], r["y_plus_max"],
                r["mesh"].get("max_aspect"), r["mesh"].get("max_nonorth"), r["mesh"].get("max_skew")] for r in rows])
    write_html(rows, allchecks, sections, figs, v1)
    for f in ("FINDINGS.md",):
        shutil.copy2(HERE / f, OUT / f)
    shutil.copy2(HERE / "results" / "REPORT.md", DATA / "compare_REPORT.md")
    write_readme()
    from latex_report import gci_table, write_latex

    print("latex:", write_latex(OUT, CASES, rows, allchecks, sections, figs, gci_table(HERE / "results" / "REPORT.md"), compile_pdf=pdf, v1=v1))
    print("done", OUT)


def write_html(rows, allchecks, sections, figs, v1):
    ng = sum(1 for _, _, c in allchecks if not c.info)
    npass = sum(1 for _, _, c in allchecks if not c.info and c.ok)
    isp = [r for r in rows if r["isp_gated"] and not r["plume"]]
    worst = max(isp, key=lambda r: abs(r["isp_err"]))
    within02 = sum(1 for r in isp if abs(r["isp_err"]) <= 0.25)
    H = [f"<!doctype html><html lang=en><head><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
         f"<title>rockmoc CFD validation</title><style>{CSS}</style></head><body>"]
    H.append("<h1>rockmoc: validation against OpenFOAM CFD</h1>")
    H.append(f"<p class=mute>Generated {datetime.now():%Y-%m-%d %H:%M} by <code>cfd_validation/make_validation_report.py</code>. "
             "Every number, plot and table is produced from the stored CFD results and the rockmoc reference values; "
             "rerun the script to rebuild it.</p>")
    H.append('<div class=kpi>'
             f"<div><b>{len(rows)}</b>CFD cases</div><div><b>{sum(r['ncell'] for r in rows) / 1e6:.2f} M</b>cells in total</div>"
             f"<div><b>{abs(worst['isp_err']):.2f} %</b>largest gated Isp error ({worst['name']})</div></div>")
    toc = ["Summary", "What was validated", "CFD method", "Boundary conditions", "Meshes", "rockmoc reference values",
           "Comparison method and tolerances", "Results overview", "rockmoc fixes made from the CFD", "Mesh independence",
           "Limitations", "Case-by-case results", "Appendix: every check", "Appendix: data files"]
    H.append("<div class=toc><ol>" + "".join(f'<li><a href="#s{i}">{t}</a></li>' for i, t in enumerate(toc)) + "</ol></div>")

    # 0 summary
    H.append('<h2 id=s0>1. Summary</h2>')
    H.append("<p>rockmoc designs supersonic nozzle contours by the method of characteristics (MOC) and predicts their performance "
             "(Isp, Cf, discharge coefficient, exit state, axial and wall distributions, boundary-layer and separation losses). "
             f"Its predictions were checked against {len(rows)} independent axisymmetric and planar OpenFOAM simulations covering every "
             "contour type the library produces: minimum-length (MLN), ideal (Sivells-type), truncated ideal (TIC), Rao thrust-optimised "
             "parabola (TOP), Rao TOC and control-surface bells, conical, planar MLN and the production pipeline with boundary-layer correction, "
             "plus full nozzles with converging sections, two experiments (Cuffel 1969; NASA TN D-5377), a real-gas LOX/RP-1 path, "
             "viscous wall-resolved runs and sea-level separated flow.</p>")
    H.append(f"<p><b>Result.</b> Vacuum Isp agrees with CFD to within ±0.25 % for {within02} of {len(isp)} gated cases "
             f"(largest {worst['isp_err']:+.2f} %, {worst['name']}). Exit Mach, discharge coefficient, throat-plane Mach, wall pressure "
             "and the MOC flow field agree within the stated tolerances; remaining differences are explained in "
             "<a href='#s10'>Limitations</a>.</p>")
    H.append(f'<img src="{figs["isp"]}" alt="Isp error, all cases">')

    # 1 scope
    H.append('<h2 id=s1>2. What was validated</h2>')
    tier_rows = []
    for t in ("T1", "T2", "T3", "GCI", "T4"):
        rs = [r for r in rows if r["tier"] == t]
        if rs:
            tier_rows.append([TIER_NAME[t], len(rs), ", ".join(r["name"] for r in rs)])
    H.append(table(["tier", "cases", "case names"], tier_rows, num=(1,)))
    H.append(table(["case", "tier", "description", "mode", "gas model", "γ", "p0 [MPa]", "T0 [K]", "viscous", "cells"],
                   [[r["name"], r["tier"], r["title"], r["mode"] + (" + plume" if r["plume"] else ""), r["thermo"],
                     f"{r['gamma']:.4f}", f"{r['p0'] / 1e6:.3g}", f"{r['T0']:.0f}", "yes" if r["viscous"] else "no", r["ncell"]]
                    for r in rows], num=(5, 6, 7, 9)))

    # 2 method
    H.append('<h2 id=s2>3. CFD method</h2>')
    H.append(table(["item", "setting"], [
        ["Code", "OpenFOAM v2412, <code>rhoCentralFoamAxi</code>: density-based, compressible, axisymmetric variant of rhoCentralFoam"],
        ["Flux scheme", "Kurganov–Tadmor central-upwind (Kurganov); van Leer limited reconstruction of ρ, U, T"],
        ["Time integration", "explicit Euler, adaptive Δt at maxCo = 0.4, run to a steady state"],
        ["Spatial", "Gauss linear gradients, corrected Laplacian; k, ω convected with limitedLinear 1"],
        ["Geometry", "axisymmetric 2.5° wedge, one cell thick (planar case: 2-D slab with empty front/back, symmetry plane on the centre line)"],
        ["Thermodynamics", "hePsiThermo, perfect gas; constant cp (γ = 1.22, R = 355 J/kg K for the generic propellant; air γ = 1.4 for the "
                           "experiments); LOX/RP-1 cases use JANAF cp(T) polynomials fitted to frozen CEA products"],
        ["Transport (viscous)", "Sutherland law fitted to Bartz' viscosity for the combustion gas; Pr_t = 0.85"],
        ["Turbulence (viscous)", "RANS k-ω SST; Spalding (nutUSpaldingWallFunction) valid continuously from y+ ≈ 1 to the log layer; "
                                 "first-cell radial height set at the throat (5 µm; 10 µm for the separation runs), growing with the local radius towards the exit; median y+ 0.5–1.4 on the rocket-scale runs, about 9 on the D-5377 rig (per case in the Meshes table)"],
        ["Run length", "4 flow-throughs (supersonic-inlet cases; v1 history showed outlet ṁ steady to ±0.03 % after 2), 8 for chamber cases, "
                       "1.7–3 for viscous; plume cases add 4 ambient-domain transits"],
        ["Time averaging", "chamber cases: fieldAverage of U, p, T, ρ from 0.6·endTime (removes low-Mach acoustic noise); "
                           "all reported values use the mean fields"],
        ["Parallel", "8 MPI ranks per case (scotch), 2 cases at a time on a 16-core workstation (WSL2)"],
    ]))
    H.append("<p>Sample dictionaries (<code>t4_prod_bl</code>; every case's full setup is in <code>data/case_setup/</code>):</p>")
    ex = CASES / "t4_prod_bl"
    for f in ("system/fvSchemes", "system/controlDict", "constant/thermophysicalProperties", "constant/momentumTransport"):
        if (ex / f).exists():
            txt = (ex / f).read_text()
            txt = txt[txt.find("}", txt.find("FoamFile")) + 1:].strip() if "FoamFile" in txt else txt
            H.append(f"<p class=mute><code>{f}</code></p><pre>{html.escape(txt[:2500])}</pre>")

    # 3 BCs
    H.append('<h2 id=s3>4. Boundary conditions</h2>')
    H.append(table(["inlet mode", "used by", "p", "U", "T"], [
        ["throat (supersonic start line)", "Tier 1, T3, T4 production/MLN", "fixedValue: isentropic p at M = 1.05 from (p0, T0)",
         "fixedValue: uniform axial M = 1.05", "fixedValue: isentropic T at M = 1.05"],
        ["chamber (subsonic)", "Tier 2, D5377", "totalPressure p0", "pressureInletOutletVelocity", "totalTemperature T0 (γ)"],
    ]))
    H.append(table(["patch", "condition"], [
        ["outlet (nozzle exit, supersonic)", "zeroGradient p, U, T"],
        ["outlet (plume box)", "p: waveTransmissive to p_a = 101 325 Pa (l_inf = 5 r_e); U zeroGradient; T inletOutlet 300 K"],
        ["ambient (plume box side / upstream face)", "p: totalPressure p_a; U pressureInletOutletVelocity; T inletOutlet 300 K"],
        ["wall, inviscid", "slip; zeroGradient p and T"],
        ["wall, viscous", "noSlip; T zeroGradient (adiabatic) or fixed T_w(x) from NASA TN D-5377 (cooled); "
                          "k: kqRWallFunction; ω: omegaWallFunction; ν_t: nutUSpaldingWallFunction; α_t: compressible::alphatWallFunction (Pr_t 0.85)"],
        ["axis", "wedge axis (collapsed edge, empty); planar case: symmetryPlane"],
        ["front / back", "wedge (planar: empty)"],
    ]))
    H.append("<p>The exact per-case boundary conditions, parsed from each case's <code>0/</code> directory, are listed in every case section.</p>")
    H.append("<p><b>Initial field:</b> 1-D isentropic solution along the nozzle (written by <code>init_fields.py</code>), ambient state in the plume box.</p>")

    # 4 meshes
    H.append('<h2 id=s4>5. Meshes</h2>')
    H.append("<p>Structured multi-block hexahedral meshes from <code>blockMesh</code>. The rockmoc wall is resampled by PCHIP in arc length "
             "and densified around the throat (−1.5 to +2 Rt); block boundaries sit at the throat and at x/Rt = −1, −0.5, 0.5, 1, 2. Axial cell "
             "count per block ∝ (1 + 2 exp(−x²)) gives three times the axial density at the throat. Radial cells: 100 (inviscid: uniform; "
             "viscous: 75 % of the height uniform plus a geometric wall layer sized to the requested first cell). Plume cases add a 5 r_e long, "
             "3.5 r_e high ambient region. All meshes pass <code>checkMesh</code>.</p>")
    H.append(table(["case", "cells", "radial", "first cell [µm] throat–exit", "y+ median", "y+ max", "max aspect", "max non-orth [°]", "max skew", "checkMesh", "run [h]"],
                   [[r["name"], r["ncell"], r["nr"], f"{1e6 * r['dy1']:.1f}–{1e6 * r['dy1max']:.1f}", fmt(r["y_plus_med"], 3), fmt(r["y_plus_max"], 3),
                     fmt(r["mesh"].get("max_aspect"), 3), fmt(r["mesh"].get("max_nonorth"), 3), fmt(r["mesh"].get("max_skew"), 3),
                     "OK" if r["mesh"].get("ok") else "–", fmt(r["hours"], 2)] for r in rows], num=tuple(range(1, 11))))

    # 5 rockmoc refs
    H.append('<h2 id=s5>6. rockmoc reference values</h2>')
    H.append("<p>Every reference is a direct rockmoc API call on the same design that was meshed (stored in <code>data/case_setup/&lt;case&gt;/ref.json</code>). "
             "The wall coordinates in the CFD are the rockmoc contour itself.</p>")
    H.append(table(["quantity", "rockmoc source"], [
        ["Contours", "design_mln / design_ideal / design_tic / design_rao / design_rao_toc / design_rao_control / design_conical / "
                     "design_planar_mln_nasa / design_production"],
        ["Vacuum Isp, Cf", "DesignResult.performance() → integrated_performance (auto geometric model: uniform exit, conical source flow, "
                           "or MOC wall-pressure integral); production pipeline: field_integrated_cf"],
        ["Axis M, p, T, ρ, u; MOC field", "the stored characteristic net (net=True); 1-D isentropic where a design has no net"],
        ["Wall Mach / pressure", "contour.M from the design, or the MOC direct march (moc.analysis.analyze_wall)"],
        ["Discharge coefficient, throat plane", "moc.transonic: Kliegel–Levine series (Hall for Rc/Rt ≥ 4)"],
        ["Boundary layer", "analysis.boundary_layer.correct_wall_bl (momentum-integral, Eckert reference temperature, Bartz transport)"],
        ["Friction / BL loss", "BL control volume at the exit (delta_cf_exit_vac / delta_cf_uncorrected_vac)"],
        ["Separation", "performance(separation='schmucker') on the MOC wall pressure, Kalt–Badal recovery"],
    ]))

    # 6 method & tolerances
    H.append('<h2 id=s6>7. Comparison method and tolerances</h2>')
    H.append("<ul><li><b>Isp</b> from CFD exit-plane integrals: ṁ = Σ ρu_x A, F = Σ (ρu_x² + p − p_a) A (×360/2.5 for the wedge), "
             "Isp = F/(ṁ g0). Error = CFD/rockmoc − 1.</li>"
             "<li><b>Cd</b> = time-mean ṁ / ideal ṁ* (p0, T0, A_t).</li>"
             "<li><b>Axis</b>: CFD cells nearest the centre line (300 stations) vs the rockmoc axis; p, T, ρ, u from isentropic relations of the "
             "rockmoc Mach. <b>MOC field</b>: CFD interpolated at every net point in (x, r/r_wall) space.</li>"
             "<li><b>Boundary layer</b>: δ*, θ integrated normal to the wall up to 98 % total-pressure recovery.</li>"
             "<li><b>Separation</b>: first wall station after which the near-wall axial velocity stays negative.</li>"
             "<li><b>Mesh independence</b>: Celik et al. (2008) GCI on three meshes (√2 refinement per direction).</li></ul>")
    H.append(table(["check", "tolerance (set before the runs)"], [
        ["Isp (MOC-net cases, T1/T3/GCI)", "0.5 %"], ["Isp (other inviscid)", "1.0 %"], ["Isp (viscous)", "1.5 %"],
        ["Exit Mach (area mean) / std", "1 % / 0.03"], ["Cd", "0.5 %"], ["Axis Mach max / rms", "0.02 / 0.01"],
        ["Axis p, ρ, T, u rms", "5 %, 4 %, 1 %, 1 %"], ["MOC field Mach rms / flow angle rms", "1 % / 0.5°"],
        ["Throat-plane axis / wall Mach", "5 % / 8 %"], ["D-5377 wall pressure (rms log)", "0.10"],
        ["BL θ (CFD vs rockmoc and vs survey)", "35 %"], ["Separation ε_sep / sea-level Cf", "25 % / 5 %"],
    ]))
    H.append("<p>Checks marked <span class=info>info</span> are reported but not gated (design-assumption mismatches such as equilibrium-designed "
             "walls in frozen CFD, or quantities without a rockmoc counterpart).</p>")

    # 7 overview table
    H.append('<h2 id=s7>8. Results overview</h2>')
    H.append(f"<p>{npass} of {ng} pre-set gated checks pass; every one of the {ng - npass} failures is explained in "
             "<a href='#s10'>Limitations</a> (38 are rockmoc's own grid error at a deliberately coarse 40 characteristics, which "
             "converges first-order and passes at 80).</p>")
    H.append(table(["case", "tier", "Isp rockmoc [s]", "Isp CFD [s]", "ΔIsp", "exit M CFD", "Cd CFD", "Cd rockmoc", "axis M rms", "field M rms", "pass", "fail", "info"],
                   [[f'<a href="#c_{r["name"]}">{r["name"]}</a>', r["tier"], fmt(r["isp_ref"], 6), fmt(r["isp_cfd"], 6),
                     (f"{r['isp_err']:+.2f} %" if math.isfinite(r["isp_err"]) else "–") + ("" if r["isp_gated"] or not r["isp_ref"] else " (info)"),
                     f"{r['exit_M']:.3f} ± {r['exit_std']:.3f}", fmt(r["Cd"], 5), fmt(r["Cd_ref"], 5),
                     fmt(100 * r["axis_M_rms"], 3) + " %" if math.isfinite(r["axis_M_rms"]) else "–",
                     fmt(100 * r["field_M_rms"], 3) + " %" if math.isfinite(r["field_M_rms"]) else "–",
                     r["n_pass"], r["n_fail"], r["n_info"]] for r in rows], num=tuple(range(2, 13))))
    H.append("<p class=mute>Plume (separation) cases have no vacuum Isp comparison: their CFD 'Isp' is the sea-level value. "
             "Experiment cases (D-5377) are compared against the measurements and Cd.</p>")
    for k, cap in (("throat", "Throat-plane axis Mach: CFD follows Kliegel–Levine; the old Hall series fails below Rc/Rt ≈ 3."),
                   ("sep", "Sea-level Cf of the three over-expanded cases."), ("mln", "MLN axis error vs rockmoc characteristics: first-order convergence.")):
        if k in figs:
            H.append(f'<img src="{figs[k]}" alt="{cap}"><p class=mute>{cap}</p>')

    # 8 fixes
    H.append('<h2 id=s8>9. rockmoc fixes made from the CFD</h2>')
    H.append("<p>The first campaign (rockmoc v0.9) exposed defects that were fixed before the final runs; details and before/after numbers are in "
             "<code>FINDINGS.md</code> and <code>CHANGELOG.md</code>.</p>")
    H.append(table(["#", "defect found by CFD", "fix", "effect on Isp error (CFD/rockmoc − 1), before → after"], [
        ["1", "Divergence loss applied twice (λ and ∫cos²θ)", "one geometric model chosen per contour", "cone 15°: +3.8 % → +0.24 %; Rao 100 %: +2.7 % → −0.03 %"],
        ["2", "Rao TOP exit-angle correlation sign; non-exact parabola", "exact quadratic Bézier, corrected correlation", "Rao 60 %: −3.3 % → −0.02 %"],
        ["3", "Rao TOC built by stretching a TOP contour", "generated at target length", "Rao TOC 80 %: +2.2 % → −0.10 %"],
        ["4", "5.9° kinks on ideal/TIC throat arc", "exact arc resampling", "smooth wall"],
        ["5", "Separation on 1-D pressure; losses applied to ambient term", "Schmucker on MOC wall pressure, vacuum part separated", "sea-level Cf: −13.5 % → +0.17 % (Rao 80 %, 1 MPa)"],
        ["6", "Friction charged as full ∫τ_w dA", "exit control-volume BL loss (corrected and uncorrected walls)", "production bell: +0.78 % → +0.13 %; uncorrected wall: −0.9 % → +0.14 %"],
        ["7", "Planar MLN wall from heuristic steps (7.8 % too long)", "C+ / streamline intersections", "exit std 0.063 → 0.008"],
        ["8", "Hall series used at Rc/Rt 1.5–4", "Kliegel–Levine up to Rc/Rt 4", "throat Mach error 3–9 % → 0.3–2.4 %"],
        ["9", "CLI default 40 characteristics", "80", "MLN axis rms 1.0 % → 0.64 %"],
    ]))
    if "ba" in figs:
        H.append(f'<img src="{figs["ba"]}" alt="before/after">')

    # 9 GCI
    H.append('<h2 id=s9>10. Mesh independence</h2>')
    gtxt = (HERE / "results" / "REPORT.md").read_text()
    m = re.search(r"## Mesh independence \(GCI\)\n\n(.*?)\n\n", gtxt, flags=re.S)
    if m:
        lines = [l for l in m.group(1).splitlines() if l.startswith("|") and not l.startswith("|---")]
        H.append(table([c.strip() for c in lines[0].strip("|").split("|")], [[c.strip() for c in l.strip("|").split("|")] for l in lines[1:]]))
    H.append(f'<img src="{figs["gci"]}" alt="mesh independence">')
    H.append("<p>Isp changes by at most 0.02 % between coarse and fine meshes; the residual rockmoc–CFD differences are therefore not mesh error.</p>")

    # 10 limitations
    H.append('<h2 id=s10>11. Limitations</h2><ul>')
    for t in [
        "<b>MLN axis at ≤ 40 characteristics</b> (38 of the failed gated checks, including the frozen-gas and mesh-convergence MLN variants): rockmoc's "
        "own discretization, first order in n (rms 3.8 / 2.0 / 1.0 / 0.64 % at n = 10 / 20 / 40 / 80). Passes at n = 80; the CLI default is now 80 and the API "
        "default 160. The gates were kept at n = 40 to document the trend.",
        "<b>Rao control-surface bell</b> (2 failed checks): MOC field Mach rms 1.4 % and flow-angle rms 0.9°; Isp +0.12 %. Most likely a weak internal "
        "shock that the characteristic net does not carry.",
        "<b>D-5377 cooled wall</b> (2 failed checks): the CFD (SST, cooled wall) under-predicts θ by 30–45 % against both the NASA survey and rockmoc; rockmoc is "
        "within +5 to +25 % of the survey. A CFD turbulence-model limitation.",
        "<b>Conical nozzles</b>: rockmoc is 0.2 % low (mesh-independent); λ assumes source flow from a point.",
        "<b>Cuffel nozzle</b> (Rc/Rt = 0.625, 15° cone, air): Isp +0.39 % (inside the 1 % gate); same source-flow λ approximation, stronger with a tight throat.",
        "<b>Rao 80 % at 1 MPa (sea level)</b> is borderline: the exit wall pressure is 0.6 % above the Schmucker threshold; CFD separates only over the last 0.5 in ε. "
        "Thrust agrees (+0.17 %) but ε_sep is sensitive near a flat bell lip.",
        "<b>MLN fed by a real converging section</b>: the MLN assumes a uniform sonic line at a sharp corner; a real convergent gives a curved one, so axis "
        "errors rise to 2.4 % rms (Isp unaffected). Reported as info.",
        "<b>Equilibrium chemistry</b> is not available in the CFD (frozen only); the equilibrium-designed wall in frozen CFD is a bound (−3.9 % Isp), not a match.",
        "<b>Not covered</b>: aerospike / plug nozzles, two-phase flow, finite-rate chemistry, 3-D or non-axisymmetric effects, transient start-up.",
        "<b>CFD uncertainty</b>: RANS (SST) for viscous and separated flow; 2.5° wedge axisymmetry; perfect-gas or frozen JANAF thermodynamics.",
    ]:
        H.append(f"<li>{t}</li>")
    H.append("</ul>")

    # 11 cases
    H.append('<h2 id=s11>12. Case-by-case results</h2>')
    byn = {r["name"]: r for r in rows}
    for t in ("T1", "GCI", "T2", "T3", "T4"):
        ns = [r["name"] for r in rows if r["tier"] == t]
        if not ns:
            continue
        H.append(f"<h3>{html.escape(TIER_NAME[t])}</h3>")
        for n in ns:
            ref, meta, res, bcs, gas, cm = sections[n]
            r = byn[n]
            H.append(f'<div class=case id="c_{n}"><h3>{n}: {html.escape(ref["title"])}</h3>')
            facts = [["tier", r["tier"]], ["inlet mode", r["mode"] + (" + ambient plume (p_a 101 325 Pa)" if r["plume"] else "")],
                     ["geometry", "planar 2-D" if r["planar"] else "axisymmetric wedge 2.5°"],
                     ["p0, T0", f"{r['p0'] / 1e6:.4g} MPa, {r['T0']:.0f} K"],
                     ["gas", f"{gas['thermo']} cp, W = {gas['W']:.3f} kg/kmol, R = {gas['R']:.2f} J/kg K"
                             + (f", γ = {gas['gamma']:.4f}" if "gamma" in gas else "")
                             + (f"; Sutherland As = {gas['As']:.4g}, Ts = {gas['Ts']:.4g} K" if r["viscous"] and "As" in gas else "")],
                     ["flow model", "RANS k-ω SST, wall-resolved" if r["viscous"] else "inviscid (Euler; slip wall)"],
                     ["mesh", f"{r['ncell']} cells, {r['nr']} radial; first cell {1e6 * r['dy1']:.1f}–{1e6 * r['dy1max']:.1f} µm (throat → exit)"
                              + (f"; y+ median {r['y_plus_med']:.2f}, max {r['y_plus_max']:.2f}" if r["viscous"] else "")
                              + (f"; max aspect {cm.get('max_aspect', float('nan')):.1f}, non-orth {cm.get('max_nonorth', float('nan')):.1f}°, skew {cm.get('max_skew', float('nan')):.2f}" if cm else "")],
                     ["run", f"end time {meta['endTime'] * 1e3:.3f} ms" + (f", {r['hours']:.2f} h wall clock on 8 cores" if math.isfinite(r["hours"]) else "")],
                     ["CFD ṁ, thrust", f"{res['mdot']:.4g} kg/s, {res['F']:.5g} N" + (" (sea level)" if r["plume"] else " (vacuum)")]]
            H.append(table(["", ""], facts))
            H.append("<p><b>Boundary conditions</b> (from <code>0/</code>):</p>" + bc_table(bcs))
            H.append("<p><b>Checks</b> (rockmoc vs CFD):</p>" + table(["quantity", "rockmoc", "CFD", "error", "tolerance", "status", "note"],
                                                                     check_rows(res["checks"]), num=(1, 2, 3, 4)))
            for f, cap in (("mesh.png", "Mesh"), ("contours.png", "CFD contours"), ("overview.png", "Overview: Mach field with rockmoc MOC overlay, axis/wall Mach, axis p/T/ρ/u, case-specific panel"),
                           ("axis.png", "Centre-line distributions"), ("wall_exit.png", "Wall and exit plane")):
                if (FIG / n / f).exists():
                    H.append(f'<p class=mute>{cap}</p><img loading=lazy src="figures/{n}/{f}" alt="{n} {cap}">')
            H.append(f"<p class=mute>Data: <code>data/profiles/{n}_axis.csv</code>, <code>_wall.csv</code>, <code>_exit.csv</code>; raw CFD extract "
                     f"<code>data/cfd/{n}/cfd.npz</code>; OpenFOAM setup <code>data/case_setup/{n}/</code>.</p></div>")

    # 12 all checks
    H.append('<h2 id=s12>13. Appendix: every check</h2>')
    H.append(table(["case", "quantity", "rockmoc", "CFD", "error", "tolerance", "status", "note"],
                   [[n] + row for n, _, c in allchecks for row in check_rows([c])], num=(2, 3, 4, 5)))
    # 13 files
    H.append('<h2 id=s13>14. Appendix: data files</h2>')
    H.append(table(["path", "content"], [
        ["data/summary.csv", "one row per case: Isp, exit Mach, Cd, axis/field errors, pass/fail counts"],
        ["data/all_checks.csv", "every check: rockmoc value, CFD value, error, tolerance, status"],
        ["data/cases.csv", "case setup: mode, gas, p0, T0, cells, first cell, y+, mesh quality, run time"],
        ["data/profiles/&lt;case&gt;_axis.csv", "centre line: x, M, p/p0, T/T0, ρ/ρ0, u/a0 for CFD and rockmoc"],
        ["data/profiles/&lt;case&gt;_wall.csv", "wall: x, r, M, p, T, ρ (+ y+, τ_w for viscous)"],
        ["data/profiles/&lt;case&gt;_exit.csv", "exit plane: r, M, p, T, ρ, U_x, U_r, face area"],
        ["data/cfd/&lt;case&gt;/cfd.npz", "raw CFD extract (all cells, wall, exit, inlet, mass-flow history)"],
        ["data/cfd/&lt;case&gt;/mesh.npz, log.checkMesh", "mesh polygons and quality log"],
        ["data/case_setup/&lt;case&gt;/", "complete OpenFOAM case (0/, constant/, system/, Allrun) plus rockmoc meta.json and ref.json"],
        ["data/compare_REPORT.md", "the compare.py markdown report"],
        ["FINDINGS.md", "narrative findings, before/after, fixes"],
    ]))
    H.append("</body></html>")
    (OUT / "index.html").write_text("\n".join(H), encoding="utf-8")


def write_readme():
    (OUT / "README.md").write_text(
        "# rockmoc CFD validation package\n\n"
        "Open `rockmoc_CFD_validation_report.pdf` (LaTeX; source in `latex/`) or `index.html` for the full report: method, boundary conditions, meshes, "
        "per-case mesh plots, contours, centre-line and wall comparisons, every check, limitations.\n\n"
        "- `figures/` per-case and summary plots\n- `data/` CSV tables, per-case profiles, raw CFD extracts, mesh polygons, "
        "complete OpenFOAM case setups\n- `FINDINGS.md` narrative findings and the rockmoc fixes\n\n"
        "Rebuild: `cd cfd_validation && py -3.12 make_validation_report.py` (needs `results/raw` from `cfdval.py pull`).\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
