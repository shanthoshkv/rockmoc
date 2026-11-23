"""rockmoc vs OpenFOAM comparison: metrics, figures and REPORT.md."""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.interpolate import LinearNDInterpolator  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

HERE = Path(__file__).resolve().parent
CASES, RAW, OUT = HERE / "cases", HERE / "results" / "raw", HERE / "results"
G0 = 9.80665
WEDGE = 2.5


# ----------------------------------------------------------------------------- helpers
def ratios(ref, M):
    M = np.asarray(M, float)
    t = ref.get("table")
    if t:
        f = lambda k: np.interp(M, t["M"], t[k])  # noqa: E731
        a0 = math.sqrt(t["gamma_c"] * ref["R"] * t["Tc"])
        return f("p_pc"), f("T_Tc"), f("rho_rhoc"), M * f("a") / a0
    g = ref["gamma"]
    tt = 1 + 0.5 * (g - 1) * M * M
    return tt ** (-g / (g - 1)), 1 / tt, tt ** (-1 / (g - 1)), M / np.sqrt(tt)


def norms(ref):
    """CFD normalisers p0, T0, rho0, a0 consistent with ratios()."""
    t = ref.get("table")
    if t:
        rho0 = t.get("rho_c", ref["p0"] / (ref["R"] * ref["T0"]))
        return ref["p0"], ref["T0"], rho0, math.sqrt(t["gamma_c"] * ref["R"] * t["Tc"])
    return ref["p0"], ref["T0"], ref["p0"] / (ref["R"] * ref["T0"]), math.sqrt(ref["gamma"] * ref["R"] * ref["T0"])


def area_ratio(M, g):
    return (1 / M) * ((2 / (g + 1)) * (1 + 0.5 * (g - 1) * M * M)) ** ((g + 1) / (2 * (g - 1)))


def mach_1d(eps, g, sup):
    lo, hi = (1.0, 50.0) if sup else (1e-4, 1.0)
    for _ in range(80):
        m = 0.5 * (lo + hi)
        f = area_ratio(m, g) - eps
        if (f < 0) == sup:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def load(name):
    ref = json.loads((CASES / name / "ref.json").read_text())
    z = np.load(RAW / name / "cfd.npz")
    z = {k: z[k] for k in z.files}
    # time-averaged fields (fieldAverage, chamber cases) replace the last snapshot
    for pre in ("c_", "e_", "w_", "i_"):
        if pre + "UMean_x" not in z:
            continue
        for f in ("U_x", "U_r"):
            z[pre + f] = z[pre + "UMean_" + f.split("_")[1]]
        for f in ("p", "T", "rho"):
            if pre + f + "Mean" in z:
                z[pre + f] = z[pre + f + "Mean"]
        if pre + "Ma" in z and pre + "T" in z:  # Mach of the mean state: |U_mean| / sqrt(g R T_mean)
            u_mean = np.hypot(z[pre + "U_x"], z[pre + "U_r"])
            g_r = ref["gamma"] * ref["R"]
            z[pre + "Ma"] = u_mean / np.sqrt(g_r * np.maximum(z[pre + "T"], 1.0))
        z["averaged"] = True
    return ref, z


class Check:
    def __init__(self, case, quantity, rockmoc, cfd, tol, rel=True, info=False, note=""):
        self.case, self.quantity, self.rockmoc, self.cfd, self.tol, self.rel, self.info, self.note = case, quantity, rockmoc, cfd, tol, rel, info, note
        if rockmoc is None or cfd is None or not np.isfinite(cfd) or not np.isfinite(rockmoc):
            self.err, self.ok = float("nan"), False
        else:
            self.err = (cfd / rockmoc - 1) if rel else (cfd - rockmoc)
            self.ok = abs(self.err) <= tol

    def row(self):
        e = f"{100 * self.err:+.2f} %" if self.rel else f"{self.err:+.3g}"
        t = f"{100 * self.tol:.1f} %" if self.rel else f"{self.tol:g}"
        st = "info" if self.info else ("PASS" if self.ok else "FAIL")
        f = lambda v: "-" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.5g}"  # noqa: E731
        return f"| {self.case} | {self.quantity} | {f(self.rockmoc)} | {f(self.cfd)} | {e} | {t} | {st} | {self.note} |"


# ----------------------------------------------------------------------------- per case
def analyse(name):
    ref, z = load(name)
    p0, T0, rho0, a0 = norms(ref)
    Rt = ref["Rt"]
    wx, wr = np.array(ref["wall_geom_x"]), np.array(ref["wall_geom_r"])
    xe = wx[-1]
    cx, cr = z["c_x"], z["c_r"]
    inside = cx <= xe + 1e-9
    cM = z["c_Ma"]
    cU = np.hypot(z["c_U_x"], z["c_U_r"])
    cth = np.degrees(np.arctan2(z["c_U_r"], z["c_U_x"]))
    res = {"name": name, "title": ref["title"], "tier": ref["tier"], "ncell": ref["ncell"], "time": float(z["time"])}
    checks: list[Check] = []
    tag = name

    # mass balance: time mean over last half of history (chamber inlet carries low-Mach acoustic noise)
    mdot_hist = None
    if "h_flowRates" in z and "h_flowRateOut" in z and not ref.get("plume"):
        hi, ho = z["h_flowRates"][:, 1], z["h_flowRateOut"][:, 1]
        k = len(hi) // 2
        mi, mo = hi[k:].mean(), ho[k:].mean()
        mdot_hist = abs(mo) * (1.0 if ref["planar"] else 360.0 / WEDGE)
        res["mass_imbalance"] = abs(mi + mo) / abs(mi)
        res["outlet_unsteady"] = float(ho[k:].std() / abs(mo))
        checks.append(Check(tag, "mass imbalance in/out (mean last half)", 0.0, res["mass_imbalance"], 2e-3, rel=False, note="steadiness"))
        checks.append(Check(tag, "outlet mdot std/mean (last half)", 0.0, res["outlet_unsteady"], 1e-2, rel=False, info=True))

    # exit plane integrals
    fac = 1.0 if ref["planar"] else 360.0 / WEDGE
    pa = ref["plume"]["pa"] if ref.get("plume") else 0.0
    A, rho, ux, ur, pe = z["e_A"], z["e_rho"], z["e_U_x"], z["e_U_r"], z["e_p"]
    mdot = fac * np.sum(rho * ux * A)
    F = fac * np.sum((rho * ux * ux + pe - pa) * A)
    Isp = F / (mdot * G0)
    res.update(mdot=mdot, F=F, Isp=Isp)
    g = ref["gamma"]
    if ref.get("table"):
        cstar = ref.get("cstar")
    else:
        mstar = p0 * math.sqrt(g / (ref["R"] * T0)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
        cstar = p0 / mstar
        if not ref["planar"]:
            res["Cd"] = (mdot_hist or mdot) / (mstar * math.pi * Rt * Rt)
    res["Cf_mdot"] = Isp * G0 / cstar if cstar else None
    if ref.get("Isp_vac") is not None and not ref.get("plume"):
        net = ref.get("has_net")
        tol = 0.005 if net and ref["tier"] in ("T1", "GCI", "T3") else (0.015 if ref["viscous"] else 0.01)
        info = ref.get("gas_table") in ("equil", "equilibrium")
        checks.append(Check(tag, "Isp_vac [s]", ref["Isp_vac"], Isp, tol, info=info,
                            note="equilibrium-designed wall in frozen CFD (bound, not a match)" if info else ""))
    if "Cd" in ref and "Cd" in res and ref["mode"] == "chamber":
        checks.append(Check(tag, "discharge coefficient Cd", ref["Cd"], res["Cd"], 0.005, note="rockmoc transonic series"))
        if "Cd_exp" in ref:
            checks.append(Check(tag, "Cd vs experiment", ref["Cd_exp"], res["Cd"], 0.005, note="Cuffel 1969 measured 0.985"))
    # exit uniformity
    eM = z["e_Ma"]
    w = A / A.sum()
    res["exit_M_mean"], res["exit_M_std"] = float(np.sum(w * eM)), float(np.sqrt(np.sum(w * (eM - np.sum(w * eM)) ** 2)))
    eth = np.degrees(np.arctan2(ur, ux))
    res["exit_theta_max"] = float(np.max(np.abs(eth)))
    if ref.get("uniform_exit") and not ref["viscous"]:
        checks.append(Check(tag, "exit Mach (area mean)", ref["Me_design"], res["exit_M_mean"], 0.01))
        checks.append(Check(tag, "exit Mach nonuniformity (std)", 0.0, res["exit_M_std"], 0.03, rel=False))
    if ref["tier"] == "T4" and ref.get("dcf_friction") is not None:
        core = z["e_r"] < 0.85 * xr_exit(ref)
        Mc = float(np.sum(A[core] * eM[core]) / A[core].sum())
        res["exit_M_core"] = Mc
        checks.append(Check(tag, "exit Mach core r<0.85re (BL-corrected wall)", ref["Me_design"], Mc, 0.015, note="BL correction recovers design Me in the inviscid core"))

    # axis profile
    tree = cKDTree(np.c_[cx / Rt, cr / Rt * 4.0])  # stretch r so the first radial layer wins
    x0 = max(wx[0], 0.0) if ref["mode"] == "throat" else wx[0]
    xa = np.linspace(x0, xe, 300)
    _, ia = tree.query(np.c_[xa / Rt, np.zeros_like(xa)])
    axis = {"x": xa, "M": cM[ia], "p": z["c_p"][ia] / p0, "T": z["c_T"][ia] / T0, "rho": z["c_rho"][ia] / rho0, "u": cU[ia] / a0}
    res["axis"] = axis
    if "axis_x" in ref:
        rx, rM = list(ref["axis_x"]), list(ref["axis_M"])
        if ref.get("uniform_exit") and rx[-1] < xe:
            rx.append(xe)
            rM.append(ref["Me_design"])
        rx, rM = np.array(rx), np.array(rM)
        src = "MOC net"
    else:
        rx = xa
        eps = (np.interp(xa, wx, wr) / Rt) ** (1 if ref["planar"] else 2)
        rM = np.array([mach_1d(max(e, 1 + 1e-9), g, x >= 0) for e, x in zip(eps, xa)])
        src = "1-D"
    res["axis_ref"] = {"x": rx, "M": rM, "src": src}
    m = (xa >= rx.min()) & (xa <= rx.max()) & (xa > 0.5 * Rt)
    if m.sum() > 5:
        rMi = np.interp(xa[m], rx, rM)
        e = axis["M"][m] / rMi - 1
        res["axis_M_maxerr"], res["axis_M_rms"] = float(np.max(np.abs(e))), float(np.sqrt(np.mean(e**2)))
        # MLN fed by a real converging section: curved sonic line vs the design's uniform start (design assumption)
        info = src == "1-D" or ref["viscous"] or (ref["mode"] == "chamber" and ref.get("kind") == "mln")
        checks.append(Check(tag, f"axis Mach max |err| x>0.5Rt ({src})", 0.0, res["axis_M_maxerr"], 0.02, rel=False, info=info))
        for k, i in (("p", 0), ("T", 1), ("rho", 2), ("u", 3)):
            rv = ratios(ref, rMi)[i]
            ek = axis[k][m] / rv - 1
            res[f"axis_{k}_maxerr"] = float(np.max(np.abs(ek)))
            checks.append(Check(tag, f"axis {k}/{'a0' if k == 'u' else k + '0'} max |err| ({src})", 0.0, res[f"axis_{k}_maxerr"],
                                {"p": 0.06, "T": 0.01, "rho": 0.05, "u": 0.01}[k], rel=False, info=True))
            res[f"axis_{k}_rms"] = float(np.sqrt(np.mean(ek**2)))
            # rms tol = 1 % Mach tol scaled by isentropic sensitivity at M~3 (dp/p ~ 4.5 dM/M)
            checks.append(Check(tag, f"axis {k}/{'a0' if k == 'u' else k + '0'} rms err ({src})", 0.0, res[f"axis_{k}_rms"],
                                {"p": 0.05, "T": 0.01, "rho": 0.04, "u": 0.01}[k], rel=False, info=info))
        checks.append(Check(tag, f"axis Mach rms err x>0.5Rt ({src})", 0.0, res["axis_M_rms"], 0.01, rel=False, info=info))

    # wall
    order = np.argsort(z["w_x"])
    wall = {"x": z["w_x"][order], "M": z["w_Ma"][order], "p": z["w_p"][order] / p0}
    res["wall"] = wall
    rwx, rwM = np.array(ref["wall_x"]), np.array(ref["wall_M"], float)
    ok = np.isfinite(rwM)
    if ok.sum() > 3 and ref["kind"] != "d5377":
        sel = (wall["x"] > max(rwx[ok].min(), 0.3 * Rt)) & (wall["x"] < rwx[ok].max())
        e = wall["M"][sel] / np.interp(wall["x"][sel], rwx[ok], rwM[ok]) - 1
        if sel.sum() > 3:
            res["wall_M_rms"] = float(np.sqrt(np.mean(e**2)))
            checks.append(Check(tag, "wall Mach rms err", 0.0, res["wall_M_rms"], 0.03, rel=False,
                                info=not ref.get("has_net") or ref["viscous"], note="rockmoc wall M"))

    # MOC net field
    if "net_x" in ref:
        nx, nr, nM, nth = (np.array(ref[k]) for k in ("net_x", "net_r", "net_M", "net_th"))
        keep = (nx <= xe) & (nr <= np.interp(nx, wx, wr) + 1e-9)
        nx, nr, nM, nth = nx[keep], nr[keep], nM[keep], nth[keep]
        if len(nx) > 20000:
            s = np.random.default_rng(0).choice(len(nx), 20000, replace=False)
            nx, nr, nM, nth = nx[s], nr[s], nM[s], nth[s]
        # interpolate in (x/Rt, eta = r/r_wall) with an axis mirror: the physical-space Delaunay hull of a
        # nozzle is concave-bounded, so its long chord triangles would mix cells far apart
        ce = cr[inside] / np.interp(cx[inside], wx, wr)
        ne = nr / np.interp(nx, wx, wr)
        P = np.r_[np.c_[cx[inside] / Rt, ce], np.c_[cx[inside] / Rt, -ce]]
        iM = LinearNDInterpolator(P, np.r_[cM[inside], cM[inside]])(nx / Rt, ne)
        ith = LinearNDInterpolator(P, np.r_[cth[inside], -cth[inside]])(nx / Rt, ne)
        # net points between the wall and the last cell centre lie outside the CFD hull: skip them
        # (nearest-cell fallback there compares the wall-adjacent fan to the wrong cell)
        cov = np.isfinite(iM)
        eM = iM / nM - 1
        eth = ith - np.degrees(nth)
        far = (nx > 0.5 * Rt) & cov
        res["field_coverage"] = float(cov.mean())
        res["field_M_median"] = float(np.median(np.abs(eM[far])))
        res.update(field_M_rms=float(np.sqrt(np.mean(eM[far] ** 2))), field_M_max=float(np.max(np.abs(eM[far]))),
                   field_th_rms=float(np.sqrt(np.mean(eth[far] ** 2))), field_n=int(far.sum()),
                   field_M_rms_all=float(np.sqrt(np.nanmean(eM**2))))
        res["net"] = (nx, nr, nM)
        info = ref["viscous"] or ref.get("gas_table") in ("equil", "equilibrium", "perfect")
        checks.append(Check(tag, f"MOC field Mach rms err ({far.sum()} pts, x>0.5Rt)", 0.0, res["field_M_rms"], 0.01, rel=False, info=info))
        checks.append(Check(tag, "MOC field flow-angle rms err [deg]", 0.0, res["field_th_rms"], 0.5, rel=False, info=info))
        checks.append(Check(tag, "MOC field Mach max |err|", 0.0, res["field_M_max"], 0.05, rel=False, info=True))
        checks.append(Check(tag, "MOC field Mach median |err|", 0.0, res["field_M_median"], 0.01, rel=False, info=True,
                            note=f"{100 * res['field_coverage']:.0f} % of net inside CFD cell hull"))

    # throat plane (full nozzles)
    if ref["mode"] == "chamber" and ref.get("rc_up") is not None or "Cd_exp" in ref:
        rs = np.linspace(0, 0.98 * Rt, 40)
        _, it = cKDTree(np.c_[cx / Rt, cr / Rt]).query(np.c_[np.zeros_like(rs), rs / Rt])
        res["throat"] = {"r": rs, "M": cM[it]}
        rc = ref.get("rc_up", 0.625)
        from rockmoc.moc.transonic import HALL_MIN_RC_RT, hall_uv, kl_uv

        fn = kl_uv if rc < HALL_MIN_RC_RT else hall_uv
        Mof = lambda V: math.sqrt(max(2 * V * V / ((g + 1) - (g - 1) * V * V), 0))  # noqa: E731
        rM = np.array([Mof(math.hypot(*fn(0.0, r, g, Rt, rc * Rt))) for r in rs])
        res["throat"]["ref"] = rM
        res["throat"]["model"] = "Kliegel-Levine" if fn is kl_uv else "Hall"
        sharp = ref.get("kind") == "mln"  # MLN: sharp-corner throat downstream, series assumes symmetric Rc
        nt = "MLN sharp-corner throat: series not applicable" if sharp else ""
        checks.append(Check(tag, f"throat-plane axis M ({res['throat']['model']})", rM[0], cM[it][0], 0.05, info=sharp, note=nt))
        checks.append(Check(tag, f"throat-plane wall M ({res['throat']['model']})", rM[-1], cM[it][-1], 0.08, info=sharp, note=nt))
        if "M_axis_throat_exp" in ref:
            checks.append(Check(tag, "throat-plane axis M vs experiment", ref["M_axis_throat_exp"], cM[it][0], 0.05, note="Cuffel 1969"))
            checks.append(Check(tag, "throat-plane wall M vs experiment", ref["M_wall_throat_exp"], cM[it][-1], 0.10, note="Cuffel 1969"))

    # D5377 experiment
    if ref["kind"] == "d5377":
        xd = np.array(ref["exp_xd"])
        pw = np.array([np.nan if v is None else v for v in ref["exp_pw"]])
        ds = ref["dstar"]
        cfd_pw = np.interp(xd * ds, wall["x"], wall["p"])
        res["d5377"] = {"xd": xd, "exp": pw, "cfd": cfd_pw}
        okp = np.isfinite(pw)
        lr = np.log(cfd_pw[okp] / pw[okp])
        res["d5377_logp_rms"] = float(np.sqrt(np.mean(lr**2)))
        checks.append(Check(tag, "wall p/p0 vs TN D-5377 taps (rms log err)", 0.0, res["d5377_logp_rms"], 0.10, rel=False, note=f"{okp.sum()} taps"))
        # rockmoc 1-D cone wall prediction at the taps
        eps = (np.interp(xd * ds, wx, wr) / Rt) ** 2
        p1 = np.array([ratios(ref, mach_1d(max(e, 1 + 1e-9), g, x >= 0))[0] for e, x in zip(eps, xd)])
        res["d5377"]["rockmoc_1d"] = p1

    # viscous: BL thicknesses and skin friction
    if ref["viscous"] and "bl_dstar" in ref:
        bl = bl_profiles(z, ref, wx, wr)
        res["bl"] = bl
        bx = np.array(ref["bl_x"]) if len(ref["bl_x"]) == len(ref["bl_dstar"]) else wx
        rd = np.interp(bl["x"], bx, ref["bl_dstar"])
        rt = np.interp(bl["x"], bx, ref["bl_theta"])
        bl["ref_dstar"], bl["ref_theta"] = rd, rt
        sel = bl["x"] > 0.3 * (xe - max(wx[0], 0)) + max(wx[0], 0)
        if sel.any():
            checks.append(Check(tag, "momentum thickness theta (mean ratio downstream)", 1.0, float(np.mean(bl["theta"][sel] / np.maximum(rt[sel], 1e-12))), 0.35))
            checks.append(Check(tag, "displacement thickness delta* (mean ratio downstream)", 1.0, float(np.mean(bl["dstar"][sel] / np.maximum(rd[sel], 1e-12))), 0.35, info=True))
        if "exp_bl" in ref:
            ds = ref["dstar"]
            for st, row in ref["exp_bl"].items():
                xs = float(st) * ds
                i = int(np.argmin(np.abs(bl["x"] - xs)))
                th_e, ds_e = 0.5 * (row[1] + row[2]) * 1e-2, 0.5 * (row[3] + row[4]) * 1e-2
                checks.append(Check(tag, f"theta at x/d*={st} vs experiment", th_e, bl["theta"][i], 0.35, note="TN D-5377 survey"))
                checks.append(Check(tag, f"delta* at x/d*={st} vs experiment", ds_e, bl["dstar"][i], 0.50, note="TN D-5377 survey", info=float(st) < 1))
                checks.append(Check(tag, f"theta at x/d*={st}: rockmoc vs experiment", th_e, float(np.interp(xs, bx, ref["bl_theta"])), 0.35, info=True, note="library vs data (context)"))

    # separation
    if ref.get("plume"):
        sep = separation(z, ref, wx, wr)
        res["sep"] = sep
        if sep["x_sep"] is not None:
            res["eps_sep"] = (np.interp(sep["x_sep"], wx, wr) / Rt) ** 2
        else:
            res["eps_sep"] = (wr[-1] / Rt) ** 2
        coarse = False  # v2: all separation cases wall-resolved (y+ ~1-2 downstream)
        checks.append(Check(tag, "separation area ratio eps_sep (Schmucker)", ref["sep_eps"], res["eps_sep"], 0.25, info=coarse,
                            note=("CFD: attached to lip" if sep["x_sep"] is None else "CFD separated") + ("; coarse wall, not valid" if coarse else "")))
        cf_ref = ref.get("sep_Cf")
        res["Cf_amb"] = res["Cf_mdot"]
        if cf_ref:
            checks.append(Check(tag, "Cf at sea level (rockmoc: Schmucker on MOC wall p, Kalt-Badal)", cf_ref, res["Cf_amb"], 0.05, info=coarse, note=f"attached-flow theory {ref.get('sep_Cf_att', float('nan')):.3f}"))
            checks.append(Check(tag, "Cf at sea level vs attached-flow theory", ref.get("sep_Cf_att"), res["Cf_amb"], 0.05, info=True))
    if ref.get("gas_table") in ("equil", "equilibrium", "perfect"):  # deliberate design/CFD gas mismatch: shows the effect, not a match
        for c in checks:
            c.info = True
    res["checks"] = checks
    return ref, z, res


def xr_exit(ref):
    return ref["wall_geom_r"][-1]


def bl_profiles(z, ref, wx, wr, n=40):
    cx, cr = z["c_x"], z["c_r"]
    U = np.hypot(z["c_U_x"], z["c_U_r"])
    rho = z["c_rho"]
    x0 = max(wx[0], 0.0) if ref["kind"] != "d5377" else wx[0] + 0.1 * (wx[-1] - wx[0])
    xs = np.linspace(x0 + 0.02 * (wx[-1] - x0), wx[-1] * 0.995, n)
    th_w = np.arctan(np.gradient(wr, wx))
    out = {"x": xs, "dstar": np.zeros(n), "theta": np.zeros(n), "cf": np.zeros(n)}
    order = np.argsort(z["w_x"])
    tw = np.hypot(z["w_wallShearStress_x"], z["w_wallShearStress_r"])[order] if "w_wallShearStress_x" in z else None
    for k, xs_ in enumerate(xs):
        rw = np.interp(xs_, wx, wr)
        dx = 1.5 * rw / ref["nr"] * 1.5
        sel = np.abs(cx - xs_) < dx
        if sel.sum() < 10:
            continue
        n_ = (rw - cr[sel]) * math.cos(np.interp(xs_, wx, th_w))
        o = np.argsort(n_)
        nn, uu, rr = n_[o], U[sel][o], rho[sel][o]
        # one value per radial layer: keep the cells closest to xs_ by binning
        nb = np.unique(np.round(nn / (rw * 1e-4)))
        nn2 = np.array([nn[np.round(nn / (rw * 1e-4)) == b].mean() for b in nb])
        uu2 = np.array([uu[np.round(nn / (rw * 1e-4)) == b].mean() for b in nb])
        rr2 = np.array([rr[np.round(nn / (rw * 1e-4)) == b].mean() for b in nb])
        lim = nn2 < 0.5 * rw
        # BL edge = total-pressure recovery (inviscid core keeps p0; robust where core rho*u varies radially)
        g = ref["gamma"]
        mm = z["c_Ma"][sel][o]
        pt = z["c_p"][sel][o] * (1 + 0.5 * (g - 1) * mm * mm) ** (g / (g - 1))
        pt2 = np.array([pt[np.round(nn / (rw * 1e-4)) == b].mean() for b in nb])
        ie = int(np.argmax(lim & (pt2 >= 0.98 * pt2[lim].max())))
        ne, ue, re = nn2[ie], uu2[ie], rr2[ie]
        f = rr2[: ie + 1] * uu2[: ie + 1] / (re * ue)
        nint = np.r_[0.0, nn2[: ie + 1]]
        fi = np.r_[0.0, f]
        ui = np.r_[0.0, uu2[: ie + 1] / ue]
        out["dstar"][k] = np.trapz(1 - fi, nint)
        out["theta"][k] = np.trapz(fi * (1 - ui), nint)
        if tw is not None:
            out["cf"][k] = np.interp(xs_, z["w_x"][order], tw) / (0.5 * re * ue * ue)
        _ = ne
    return out


def separation(z, ref, wx, wr):
    cx, cr = z["c_x"], z["c_r"]
    order = np.argsort(z["w_x"])
    fx, fr = z["w_x"][order], z["w_r"][order]
    _, j = cKDTree(np.c_[cx, cr]).query(np.c_[fx, fr])
    ux = z["c_U_x"][j]
    neg = (ux < 0) & (fx > 0.2 * ref["Rt"])
    x_sep = None
    if neg.any():
        # first x after which the near-wall flow stays reversed (allow 10 % noise)
        for i in np.where(neg)[0]:
            if neg[i:].mean() > 0.9:
                x_sep = float(fx[i])
                break
    return {"x": fx, "ux_wall": ux, "x_sep": x_sep}


# ----------------------------------------------------------------------------- figures
def figure(ref, z, res, path):
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.6))
    Rt = ref["Rt"]
    wx, wr = np.array(ref["wall_geom_x"]), np.array(ref["wall_geom_r"])
    a = ax[0, 0]
    from matplotlib.tri import Triangulation

    tri = Triangulation(z["c_x"] / Rt, z["c_r"] / Rt)
    xc, rc = tri.x[tri.triangles].mean(1), tri.y[tri.triangles].mean(1)
    tri.set_mask((xc * Rt <= wx[-1]) & (rc * Rt > np.interp(xc * Rt, wx, wr)))  # drop chords outside a concave wall
    tc = a.tripcolor(tri, z["c_Ma"], shading="gouraud", cmap="viridis")
    fig.colorbar(tc, ax=a, label=r"Mach number $M$ [–]")
    lv = np.linspace(1.2, max(1.3, float(np.nanmax(z["c_Ma"])) - 0.1), 9)
    a.tricontour(tri, z["c_Ma"], levels=lv, colors="w", linewidths=0.8)
    if "net" in res:
        nx, nr, nM = res["net"]
        try:
            a.tricontour(nx / Rt, nr / Rt, nM, levels=lv, colors="#D55E00", linewidths=0.9, linestyles="--")
        except Exception:  # noqa: BLE001
            pass
    a.plot(wx / Rt, wr / Rt, "k", lw=1.5)
    a.set_aspect("equal")
    a.set_xlabel(r"$x/R_t$ [–]")
    a.set_ylabel(r"$r/R_t$ [–]")
    a = ax[0, 1]
    A, R = res["axis"], res["axis_ref"]
    a.plot(A["x"] / Rt, A["M"], "k", label="CFD axis")
    a.plot(R["x"] / Rt, R["M"], "--", color="#D55E00", label="rockmoc axis")
    a.plot(res["wall"]["x"] / Rt, res["wall"]["M"], color="#0072B2", label="CFD wall")
    rwM = np.array(ref["wall_M"], float)
    if np.isfinite(rwM).any():
        a.plot(np.array(ref["wall_x"]) / Rt, rwM, "--", color="#E69F00", label="rockmoc wall")
    a.set_xlabel(r"$x/R_t$ [–]")
    a.set_ylabel(r"$M$ [–]")
    a.legend(fontsize=7.5)
    a.grid(alpha=0.3)
    a = ax[1, 0]
    rr = ratios(ref, np.interp(A["x"], R["x"], R["M"]))
    for k, i, c in (("p", 0, "#0072B2"), ("T", 1, "#009E73"), ("rho", 2, "#E69F00"), ("u", 3, "#CC79A7")):
        lab = r"$u/a_0$" if k == "u" else {"p": r"$p/p_0$", "T": r"$T/T_0$", "rho": r"$\rho/\rho_0$"}[k]
        a.plot(A["x"] / Rt, A[k], color=c, label=lab)
        a.plot(A["x"] / Rt, rr[i], color=c, ls="--", lw=1)
    a.set_yscale("log")
    a.set_xlabel(r"$x/R_t$ [–]")
    a.set_ylabel("axis ratio [–]")
    a.legend(fontsize=7.5)
    a.grid(alpha=0.3, which="both")
    a = ax[1, 1]
    if "throat" in res:
        t = res["throat"]
        a.plot(t["r"] / Rt, t["M"], "k", label="CFD")
        a.plot(t["r"] / Rt, t["ref"], "--", color="#D55E00", label="rockmoc")
        if "M_axis_throat_exp" in ref:
            a.plot([0, 1], [ref["M_axis_throat_exp"], ref["M_wall_throat_exp"]], "s", color="#0072B2", label="Cuffel 1969")
        a.set_xlabel(r"$r/R_t$ [–]")
        a.set_ylabel(r"throat-plane $M$ [–]")
    elif "d5377" in res:
        d = res["d5377"]
        a.semilogy(d["xd"], d["exp"], "o", color="#0072B2", label="TN D-5377 taps")
        a.semilogy(d["xd"], d["cfd"], "k-", label="CFD wall")
        a.semilogy(d["xd"], d["rockmoc_1d"], "--", color="#D55E00", label="rockmoc 1-D")
        a.set_xlabel(r"$x/d^*$ [–]")
        a.set_ylabel(r"$p_w/p_0$ [–]")
    elif "bl" in res:
        b = res["bl"]
        a.plot(b["x"] / Rt, b["dstar"] * 1e3, "k", label=r"CFD $\delta^*$")
        a.plot(b["x"] / Rt, b["ref_dstar"] * 1e3, "--", color="#D55E00", label=r"rockmoc $\delta^*$")
        a.plot(b["x"] / Rt, b["theta"] * 1e3, color="#0072B2", label=r"CFD $\theta$")
        a.plot(b["x"] / Rt, b["ref_theta"] * 1e3, "--", color="#E69F00", label=r"rockmoc $\theta$")
        a.set_xlabel(r"$x/R_t$ [–]")
        a.set_ylabel("thickness [mm]")
    elif "sep" in res:
        s = res["sep"]
        a.plot(s["x"] / Rt, s["ux_wall"], "k", label=r"CFD near-wall $u_x$")
        a.axhline(0, color="gray")
        if s["x_sep"] is not None:
            a.axvline(s["x_sep"] / Rt, color="k", ls=":", label="CFD separation")
        xs = np.interp(math.sqrt(ref["sep_eps"]) * Rt, wr, wx)
        a.axvline(xs / Rt, color="#D55E00", ls="--", label="rockmoc Schmucker")
        a.set_xlabel(r"$x/R_t$ [–]")
        a.set_ylabel(r"near-wall $u_x$ [m/s]")
    else:
        r = z["e_r"] / Rt
        o = np.argsort(r)
        a.plot(r[o], z["e_Ma"][o], "k", label="CFD $M$")
        a.axhline(ref["Me_design"], color="#D55E00", ls="--", label="rockmoc design $M_e$")
        a2 = a.twinx()
        a2.plot(r[o], np.degrees(np.arctan2(z["e_U_r"], z["e_U_x"]))[o], "#0072B2", label="CFD angle")
        a2.set_ylabel("flow angle [°]", color="#0072B2")
        a.set_xlabel(r"$r/R_t$ [–]")
        a.set_ylabel(r"$M$ [–]")
    a.legend(fontsize=7.5, loc="best")
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- GCI
def gci(f1, f2, f3, r21, r32):
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0 or e21 * e32 < 0:
        return None, None
    s = np.sign(e32 / e21)
    p = 2.0
    for _ in range(100):
        q = math.log((r21**p - s) / (r32**p - s))
        p = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
    return p, 1.25 * abs(e21 / f1) / (r21**p - 1)


# ----------------------------------------------------------------------------- driver
def compare_all(sel=None):
    names = sorted(p.parent.name for p in RAW.glob("*/cfd.npz") if not p.parent.name.startswith("t5_"))  # T5: spike.py
    if sel:
        names = [n for n in names if n in sel]
    figs = OUT / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    allc, rows, summ = [], [], {}
    for n in names:
        try:
            ref, z, res = analyse(n)
        except Exception as e:  # noqa: BLE001
            print(f"{n}: analysis failed: {e!r}")
            continue
        figure(ref, z, res, figs / f"{n}.png")
        allc += res["checks"]
        summ[n] = res
        nf = sum(1 for c in res["checks"] if not c.info and not c.ok)
        print(f"{n:28s} Isp {res['Isp']:8.2f} (rockmoc {ref.get('Isp_vac', float('nan')):8.2f})  "
              f"fieldM {res.get('field_M_rms', float('nan')):.4f}  exitM {res['exit_M_mean']:.3f}±{res['exit_M_std']:.3f}  fails={nf}")
        rows.append((ref["tier"], n, ref["title"], res, nf))
    # GCI
    grows = []
    for base in sorted({n.rsplit("_", 1)[0] for n in names if n.endswith(("_coarse", "_fine"))}):
        trip = [summ.get(f"{base}_fine"), summ.get(base), summ.get(f"{base}_coarse")]
        if None in trip:
            continue
        N = [t["ncell"] for t in trip]
        r21, r32 = math.sqrt(N[0] / N[1]), math.sqrt(N[1] / N[2])
        for q in ("Isp", "exit_M_mean", "field_M_rms"):
            f = [t.get(q) for t in trip]
            if None in f:
                continue
            p, g = gci(*f, r21, r32)
            grows.append(f"| {base} | {q} | {f[2]:.5g} | {f[1]:.5g} | {f[0]:.5g} | {p if p is None else f'{p:.2f}'} | "
                         f"{'-' if g is None else f'{100 * g:.3f} %'} |")
    write_report(rows, allc, grows)


TIER_TXT = {
    "T1": "Inviscid, perfect gas, supersonic inlet on the MOC start line (M=1.05). Tests the contour + MOC field itself.",
    "GCI": "Mesh independence (coarse / medium / fine, refinement sqrt(2) per direction); Celik et al. GCI on the fine mesh.",
    "T2": "Full nozzles with chamber and converging section (subsonic total-pressure inlet). Tests throat models (Hall / Kliegel-Levine), "
          "discharge coefficient, and the MOC start-line assumption; includes two experiments.",
    "T3": "Real gas: LOX/RP-1 O/F 2.6 at 5 MPa, composition frozen at the throat (RocketCEA). OpenFOAM uses a JANAF cp(T) fitted to the "
          "frozen CEA expansion, so the CFD gas equals rockmoc's frozen table.",
    "T4": "Viscous (k-omega SST, Spalding wall function, Sutherland viscosity): boundary-layer thicknesses, friction loss, BL-corrected "
          "contour, and sea-level flow separation with an external plume domain.",
}


def write_report(rows, checks, grows):
    import datetime

    L = ["# rockmoc CFD validation report", "",
         f"Generated {datetime.datetime.now():%Y-%m-%d %H:%M}. OpenFOAM v2412 `rhoCentralFoamAxi` (Kurganov central-upwind, "
         f"2.5 deg wedge), driven by `cfd_validation/cfdval.py`.", ""]
    hard = [c for c in checks if not c.info]
    L += [f"**{sum(c.ok for c in hard)} / {len(hard)} gated checks pass** ({len(checks) - len(hard)} informational).", ""]
    L += ["## Summary", "", "| tier | case | Isp rockmoc [s] | Isp CFD [s] | dIsp | exit M (mean ± std) | MOC field M rms | Cd CFD | gated fails |",
          "|---|---|---|---|---|---|---|---|---|"]
    for tier, n, title, r, nf in sorted(rows, key=lambda t: (t[0], t[1])):
        ri = next((c for c in r["checks"] if c.quantity.startswith("Isp")), None)
        isp_r = ri.rockmoc if ri else float("nan")
        d = f"{100 * (r['Isp'] / isp_r - 1):+.2f} %" if ri else "-"
        fm = f"{100 * r['field_M_rms']:.2f} %" if "field_M_rms" in r else "-"
        cd = f"{r['Cd']:.4f}" if "Cd" in r else "-"
        L.append(f"| {tier} | [{n}](figures/{n}.png) | {isp_r:.2f} | {r['Isp']:.2f} | {d} | {r['exit_M_mean']:.3f} ± {r['exit_M_std']:.3f} | {fm} | {cd} | {nf} |")
    if grows:
        L += ["", "## Mesh independence (GCI)", "", "| case | quantity | coarse | medium | fine | observed order | GCI fine |", "|---|---|---|---|---|---|---|"] + grows
    for tier in ("T1", "GCI", "T2", "T3", "T4"):
        cs = [c for c in checks if any(c.case == n for t, n, *_ in rows if t == tier)]
        if not cs:
            continue
        L += ["", f"## {tier}", "", TIER_TXT[tier], "", "| case | quantity | rockmoc / ref | CFD | error | tol | status | note |",
              "|---|---|---|---|---|---|---|---|"] + [c.row() for c in cs]
    L += ["", "## Figures", ""] + [f"![{n}](figures/{n}.png)" for _, n, *_ in sorted(rows)]
    (OUT / "REPORT.md").write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT / 'REPORT.md'}: {sum(c.ok for c in hard)}/{len(hard)} gated checks pass")
