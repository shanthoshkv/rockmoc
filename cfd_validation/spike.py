"""T5: annular aerospike (Lee / Vernacchia spike) vs OpenFOAM, inviscid, 2.5 deg wedge.

Domain (5 blocks): spike region bounded by the sonic throat line (inlet), the spike wall,
the design jet boundary r = Re and the end-of-spike station; a base block, a near-wake
block and two ambient blocks above r = Re. At the design pressure ratio the free jet
boundary from the lip is straight (p = pe, axial flow), so r = Re is the natural block line.

    py -3.12 cfd_validation/cfdval.py build t5_spike_design ...   # cases come from build_cases()
    py -3.12 cfd_validation/spike.py compare                       # -> results/spike/REPORT.md + figures
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "spike"
RAW = HERE / "results" / "raw"

# design point (ε=8 keeps the throat gap at ~6 % of Re, resolvable on the wedge)
P = dict(Pc=5e6, Tc=3200.0, gamma=1.25, molar_mass=22.0, Re=0.05, expansion_ratio=8.0)
R_CUT = 0.04  # full spike ends at R = 0.04 Re (tip region carries ~0 thrust, avoids a degenerate axis cell)


def spike_cases():
    import cfdval as cv
    from rockmoc import design_spike, design_spike_truncated

    gas = cv.Gas(P["gamma"], cv.RU / P["molar_mass"])
    full = design_spike(**P, Pa=0.0, n_points=4000)
    pe = full.Pe
    out = []
    specs = [("t5_spike_design", "Aerospike full spike, design NPR (pa = pe)", None, 1.0),
             ("t5_spike_pa3", "Aerospike full spike, overexpanded (pa = 3 pe)", None, 3.0),
             ("t5_spike_pa0.3", "Aerospike full spike, underexpanded (pa = 0.3 pe)", None, 0.3),
             ("t5_spike_trunc25", "Aerospike 25 % truncated spike, design NPR: base pressure", 0.25, 1.0)]
    for name, title, frac, kpa in specs:
        pa = kpa * pe
        res = design_spike_truncated(full, length_fraction=frac) if frac else full
        keep = res.R >= (R_CUT * P["Re"] if not frac else 0.0)
        X, R = res.X[keep], res.R[keep]
        ref_pa = design_spike(**P, Pa=pa, n_points=4000)  # rockmoc's ambient performance (design wall pressure)
        if frac:
            ref_pa = design_spike_truncated(ref_pa, length_fraction=frac)
        ref = {"spike": True, "Re": P["Re"], "pa": pa, "pe": pe, "Me_design": full.exit_mach, "delta": full.delta_rad,
               "ht": full.ht, "At": full.At, "Cf_rockmoc": ref_pa.Cf, "Isp_rockmoc": ref_pa.Isp, "mdot": full.m_dot,
               "X": X.tolist(), "Rs": R.tolist(), "M": res.M[keep].tolist(), "P_pc": (res.P[keep] / P["Pc"]).tolist(),
               "truncated": bool(frac), "kind": "aerospike", "uniform_exit": not frac}
        Mi = 1.01
        c = cv.Case(name, "T5", title, X.copy(), R.copy(), full.ht, gas, P["Pc"], P["Tc"],
                    inlet=cv.isen_inlet(gas, P["Pc"], P["Tc"], Mi), plume={"pa": pa, "Ta": 300.0, "p_init": min(pa, pe)}, nr=60, scale=0.5, ref=ref)
        c.n_ft, c.mesher, c.inlet_angle = 6.0, mesh_spike, -(0.5 * math.pi - full.delta_rad)  # throat flow turned ν_e toward axis
        out.append(c)
    return out


def _hex(b, c, n, grading, patches, bottom=None):
    """Hex from 4 (x, r) corners c = [p00, p10, p11, p01] (counter-clockwise from lower-left)."""
    hx = [b.v(*c[0], 0), b.v(*c[1], 0), b.v(*c[2], 0), b.v(*c[3], 0),
          b.v(*c[0], 1), b.v(*c[1], 1), b.v(*c[2], 1), b.v(*c[3], 1)]
    b.blocks.append(f"hex ({' '.join(map(str, hx))}) ({n[0]} {n[1]} 1) simpleGrading ({grading[0]:.6g} {grading[1]:.6g} 1)")
    if bottom is not None and len(bottom):
        for side, (i, j) in ((0, (0, 1)), (1, (4, 5))):
            pts = " ".join(f"({p[0]:.9g} {p[1]:.9g} {p[2]:.9g})" for p in (b.pt(x, r, side) for x, r in bottom))
            b.edges.append(f"polyLine {hx[i]} {hx[j]} ({pts})")
    fdef = {"xmin": (0, 4, 7, 3), "xmax": (1, 2, 6, 5), "rmin": (0, 1, 5, 4), "rmax": (3, 7, 6, 2),
            "back": (0, 3, 2, 1), "front": (4, 5, 6, 7)}
    for side, patch in {**patches, "back": "back", "front": "front"}.items():
        b.faces.setdefault(patch, []).append("(" + " ".join(str(hx[i]) for i in fdef[side]) + ")")


def mesh_spike(case, path: Path):
    import cfdval as cv

    b = cv.Block(False)
    x, r = case.wall_x, case.wall_r
    Re, s = case.ref["Re"], case.scale
    xc, rc = float(x[-1]), float(r[-1])
    nr = max(12, round(case.nr * s))
    nx1 = max(40, round(5 * case.nr * s))
    nb = max(4, round(nr * rc / Re))  # base block, ~same spacing as the right edge of block 1
    ld, rf = 2.0 * Re, 2.5 * Re
    nx2, no = max(10, round(1.2 * case.nr * s)), max(10, round(0.5 * case.nr * s))
    gx1 = 10.0  # last/first: cluster at the lip / throat where the centred fan sits
    pts = list(zip(x[1:-1:max(1, len(x) // 400)], r[1:-1:max(1, len(x) // 400)]))
    _hex(b, [(x[0], r[0]), (xc, rc), (xc, Re), (0.0, Re)], (nx1, nr), (gx1, 1.0),
         {"xmin": "inlet", "rmin": "wall"}, bottom=pts)
    _hex(b, [(xc, rc), (xc + ld, rc), (xc + ld, Re), (xc, Re)], (nx2, nr), (3.0, 1.0), {"xmax": "outlet"})
    if rc > 1e-9:
        _hex(b, [(xc, 0.0), (xc + ld, 0.0), (xc + ld, rc), (xc, rc)], (nx2, nb), (3.0, 1.0),
             {"xmin": "wall", "xmax": "outlet", "rmin": "axis"})
    _hex(b, [(0.0, Re), (xc, Re), (xc, rf), (0.0, rf)], (nx1, no), (gx1, 4.0), {"xmin": "ambient", "rmax": "ambient"})
    _hex(b, [(xc, Re), (xc + ld, Re), (xc + ld, rf), (xc, rf)], (nx2, no), (3.0, 4.0), {"xmax": "outlet", "rmax": "ambient"})
    ptypes = {"inlet": "patch", "outlet": "patch", "wall": "wall", "ambient": "patch", "axis": "empty",
              "front": "wedge", "back": "wedge"}
    b.write(path / "system" / "blockMeshDict", ptypes)
    return nx1 * (nr + no) + nx2 * (nr + nb + no), nr, 1.0


# ----------------------------------------------------------------------------- compare
def analyse(name):
    ref = json.loads((HERE / "cases" / name / "ref.json").read_text())
    z = dict(np.load(RAW / name / "cfd.npz"))
    g, R, p0, pa = ref["gamma"], ref["R"], ref["p0"], ref["pa"]
    # thrust = throat momentum + pressure flux + ∫(p - pa) dA_x on spike and base (Lee control volume)
    U = np.hypot(z["i_U_x"], z["i_U_r"])
    rho = z["i_p"] / (R * z["i_T"])
    F_in = np.sum((rho * U * z["i_U_x"] + (z["i_p"] - pa) * np.abs(z["i_U_x"]) / U) * z["i_A"])
    F_w = np.sum((z["w_p"] - pa) * np.abs(z["w_nx"]) * z["w_A"])
    A_t = np.sum(z["i_A"])
    Cf = (F_in + F_w) / (p0 * A_t)
    spike = np.abs(z["w_nx"]) < 0.999
    o = np.argsort(z["w_x"][spike])
    wx, wp = z["w_x"][spike][o], z["w_p"][spike][o] / p0
    wM = np.hypot(z["w_U_x"], z["w_U_r"])[spike][o] / np.sqrt(g * R * z["w_T"][spike][o])
    X, Pd, Md = np.array(ref["X"]), np.array(ref["P_pc"]), np.array(ref["M"])
    inside = (wx > X[0] + 0.02 * (X[-1] - X[0])) & (wx < X[-1])
    dp = np.abs(wp[inside] - np.interp(wx[inside], X, Pd)) / np.interp(wx[inside], X, Pd)
    dM = np.abs(wM[inside] - np.interp(wx[inside], X, Md)) / np.interp(wx[inside], X, Md)
    e = z["e_r"] < ref["Re"]
    eM = np.hypot(z["e_U_x"][e], z["e_U_r"][e]) / np.sqrt(g * R * z["e_T"][e])
    base = np.abs(z["w_nx"]) >= 0.999
    res = {"Cf_cfd": float(Cf), "Cf_rockmoc": ref["Cf_rockmoc"], "dCf_pct": 100 * (Cf / ref["Cf_rockmoc"] - 1),
           "wall_p_err_med_pct": 100 * float(np.median(dp)), "wall_p_err_max_pct": 100 * float(np.max(dp)),
           "wall_M_err_med_pct": 100 * float(np.median(dM)), "exit_M_mean": float(np.mean(eM)), "exit_M_std": float(np.std(eM)),
           "Me_design": ref["Me_design"], "pa_pe": ref["pa"] / ref["pe"],
           "mdot_out_in": float(abs(z["h_flowRateOut"][-1, 1] / z["h_flowRates"][-1, 1])),
           "p_base_pa": float(np.sum(z["w_p"][base] * z["w_A"][base]) / np.sum(z["w_A"][base]) / pa) if base.any() else None}
    return ref, z, res, (wx, wp, wM, X, Pd, Md)


def compare(sel=None):
    import matplotlib.pyplot as plt

    OUT.mkdir(parents=True, exist_ok=True)
    names = sorted(p.parent.name for p in RAW.glob("t5_spike*/cfd.npz") if not sel or p.parent.name in sel)
    rows = []
    for n in names:
        ref, z, res, (wx, wp, wM, X, Pd, Md) = analyse(n)
        rows.append((n, ref["title"], res))
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
        ax[0].plot(X / ref["Re"], Pd, "k-", lw=1, label="rockmoc (Lee)")
        ax[0].plot(wx / ref["Re"], wp, ".", ms=2, color="#0072B2", label="OpenFOAM")
        ax[0].axhline(ref["pa"] / ref["p0"], color="#D55E00", ls="--", lw=0.8, label="$p_a/p_c$")
        ax[0].set(yscale="log", xlabel="x / Re", ylabel="spike wall p / pc")
        ax[1].plot(X / ref["Re"], Md, "k-", lw=1)
        ax[1].plot(wx / ref["Re"], wM, ".", ms=2, color="#0072B2")
        ax[1].set(xlabel="x / Re", ylabel="spike wall Mach")
        ax[0].legend(frameon=False, fontsize=8)
        fig.suptitle(ref["title"], fontsize=10)
        fig.tight_layout()
        fig.savefig(OUT / f"{n}.png", dpi=150)
        plt.close(fig)
    lines = ["# T5 aerospike: rockmoc (Lee/Vernacchia) vs OpenFOAM rhoCentralFoamAxi (inviscid)", "",
             "| case | pa/pe | Cf CFD | Cf rockmoc | ΔCf % | wall p err med/max % | wall M err med % | exit M (design Me) | p_base/pa | ṁ out/in |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for n, _, r in rows:
        pb = f"{r['p_base_pa']:.2f}" if r["p_base_pa"] is not None else "—"
        lines.append(f"| {n} | {r['pa_pe']:.2g} | {r['Cf_cfd']:.4f} | {r['Cf_rockmoc']:.4f} | {r['dCf_pct']:+.2f} | "
                     f"{r['wall_p_err_med_pct']:.2f} / {r['wall_p_err_max_pct']:.1f} | {r['wall_M_err_med_pct']:.2f} | "
                     f"{r['exit_M_mean']:.3f}±{r['exit_M_std']:.3f} ({r['Me_design']:.3f}) | {pb} | {r['mdot_out_in']:.3f} |")
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (OUT / "results.json").write_text(json.dumps({n: r for n, _, r in rows}, indent=1))
    print("\n".join(lines))


if __name__ == "__main__":
    import sys

    if sys.argv[1:2] == ["compare"]:
        compare(sys.argv[2:])
    else:
        sys.exit(__doc__)
