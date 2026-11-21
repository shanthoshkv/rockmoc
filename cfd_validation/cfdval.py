"""rockmoc CFD validation driver (OpenFOAM v2412 on the remote PC).

    py -3.12 cfd_validation/cfdval.py build   [case ...]   # write OpenFOAM cases + rockmoc reference
    py -3.12 cfd_validation/cfdval.py launch  [case ...]   # push + start 2 lanes x 8 cores (detached)
    py -3.12 cfd_validation/cfdval.py status
    py -3.12 cfd_validation/cfdval.py pull                 # fetch cfd.npz + logs of finished cases
    py -3.12 cfd_validation/cfdval.py compare [case ...]   # rockmoc vs CFD -> results/, REPORT.md

Case = rockmoc design -> wedge (or planar) blockMesh of the wall -> rhoCentralFoamAxi
(inviscid or kOmegaSST) -> extract.py on the PC -> compare here against the rockmoc
MOC net, axis, wall, exit plane, performance, throat, BL and separation predictions.
"""
from __future__ import annotations

import copy
import json
import math
import shutil
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))
sys.path.insert(0, r"C:\Users\kvsha\Documents\Learning\openfoam scripting website\scripts")
warnings.filterwarnings("ignore")

CASES_DIR = HERE / "cases"
RES_DIR = HERE / "results"
REMOTE = "~/foamstudio/nozzles/cfdval"
G0 = 9.80665
RU = 8314.462618
WEDGE = 2.5  # deg, total wedge angle
LANES = {"A": "0,2,4,6,8,10,12,14", "B": "16,18,20,22,24,26,28,30"}  # physical cores (SMT siblings are n, n+1)
NP = 8


# ----------------------------------------------------------------------------- gas
@dataclass
class Gas:
    gamma: float
    R: float
    kind: str = "perfect"  # perfect | janaf
    janaf: list[float] | None = None  # 7 coeffs, cp/R polynomial (both ranges)
    sutherland: tuple[float, float] | None = None  # (As, Ts) for viscous runs
    table: dict | None = None  # real-gas reference arrays: M, p_pc, T_Tc, rho_rhoc, a

    def ratios(self, M):
        """p/p0, T/T0, rho/rho0, u/a0 at Mach M (isentropic, or CEA table)."""
        M = np.asarray(M, float)
        if self.table is not None:
            t = self.table
            f = lambda k: np.interp(M, t["M"], t[k])
            a0 = math.sqrt(t["gamma_c"] * self.R * t["Tc"])
            return f("p_pc"), f("T_Tc"), f("rho_rhoc"), M * f("a") / a0
        g = self.gamma
        tt = 1 + 0.5 * (g - 1) * M * M
        return tt ** (-g / (g - 1)), 1 / tt, tt ** (-1 / (g - 1)), M / np.sqrt(tt)


def ideal_isp(g, R, T0, Me):
    tt = 1 + 0.5 * (g - 1) * Me**2
    Te = T0 / tt
    ue = Me * math.sqrt(g * R * Te)
    return (ue + R * Te / ue) / G0  # uniform parallel exit, vacuum: (u + p/(rho u))


def mdot_ideal_per_area(g, R, p0, T0):
    return p0 * math.sqrt(g / (R * T0)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))


def sutherland_fit(mu_fn, T1=800.0, T2=3000.0):
    m1, m2 = float(mu_fn(T1)), float(mu_fn(T2))
    q = m1 / m2 * (T2 / T1) ** 1.5
    Ts = (T2 - q * T1) / (q - 1)
    return m1 * (T1 + Ts) / T1**1.5, Ts


# ----------------------------------------------------------------------------- case
@dataclass
class Case:
    name: str
    tier: str
    title: str
    wall_x: np.ndarray
    wall_r: np.ndarray
    Rt: float
    gas: Gas
    p0: float
    T0: float
    mode: str = "throat"  # throat (supersonic inlet at M_t) | chamber (subsonic total p/T inlet)
    inlet: tuple[float, float, float] | None = None  # (M, p, T) for throat mode
    planar: bool = False
    viscous: dict | None = None  # {"Tw": None | float | ([x],[Tw]), "dy_wall": m}
    plume: dict | None = None  # {"pa": Pa, "Ta": K}
    nr: int = 100
    ar: float = 1.5
    scale: float = 1.0
    n_ft: float = 4.0  # supersonic-inlet cases settle within ~2 flow-throughs (v1 history)
    breaks: list[float] = field(default_factory=list)
    ref: dict = field(default_factory=dict)


def isen_inlet(gas: Gas, p0, T0, M=1.05):
    if gas.table is not None:
        t = gas.table
        return (M, float(np.interp(M, t["M"], t["p_pc"])) * p0, float(np.interp(M, t["M"], t["T_Tc"])) * T0)
    pr, Tr, _, _ = gas.ratios(M)
    return (M, float(pr) * p0, float(Tr) * T0)


def converging(Rt, rc_up=1.5, half_deg=30.0, cr=4.0, lch=1.5):
    """Chamber cylinder + cone + upstream throat arc, ending at (0, Rt) with zero slope."""
    b = math.radians(half_deg)
    Rc = rc_up * Rt
    Rch = Rt * math.sqrt(cr)
    xa = -Rc * math.sin(b)
    ra = Rt + Rc * (1 - math.cos(b))
    xcone = xa - (Rch - ra) / math.tan(b)
    x0 = xcone - lch * Rt
    xs = [np.linspace(x0, xcone, 12)[:-1], np.linspace(xcone, xa, 30)[:-1], np.linspace(xa, 0.0, 40)]
    x = np.concatenate(xs)
    r = np.where(x <= xcone, Rch, np.where(x <= xa, ra + (xa - x) * math.tan(b), Rt + Rc - np.sqrt(np.maximum(Rc**2 - x**2, 0))))
    return x, r, [xcone, xa]


def join(xc, rc, xd, rd):
    return np.concatenate([xc[:-1], xd]), np.concatenate([rc[:-1], rd])


def net_ref(res, xe):
    n = res.net
    if n is None:
        return {}
    P = [p for p in n.points if p.x <= xe + 1e-9]
    out = {"net_x": [p.x for p in P], "net_r": [p.r for p in P], "net_M": [p.M for p in P], "net_th": [p.theta for p in P]}
    ax = sorted((p.x, p.M) for p in n.axis)
    if ax:
        out["axis_x"], out["axis_M"] = [a[0] for a in ax], [a[1] for a in ax]
    return out


def from_design(name, tier, title, res, gas: Gas, p0=5e6, T0=3200.0, **kw) -> Case:
    c = res.contour
    ref = {"wall_x": c.x.tolist(), "wall_M": (c.M if c.M is not None else np.full_like(c.x, np.nan)).tolist(),
           "Me_design": float(res.exit_mach_design), "kind": res.kind, "has_net": res.net is not None}
    ref.update(net_ref(res, c.x[-1]))
    try:
        perf = res.performance(pc=p0, Tc=T0)
        ref.update(Isp_vac=perf["Isp_vac_s"], Cf_vac=perf["Cf_vac"], cstar=perf["cstar_m_s"])
    except Exception as e:  # noqa: BLE001
        ref["perf_error"] = str(e)
    ref["uniform_exit"] = res.kind in ("mln", "ideal")
    case = Case(name, tier, title, c.x.copy(), c.r.copy(), res.spec.throat_radius, gas, p0, T0, ref=ref, **kw)
    if case.mode == "throat" and case.inlet is None:
        case.inlet = isen_inlet(gas, p0, T0, res.spec.throat_mach)
    return case


def full_nozzle(case: Case, rc_up=1.5, half_deg=30.0, cr=4.0) -> Case:
    xc, rc, br = converging(case.Rt, rc_up, half_deg, cr)
    case.wall_x, case.wall_r = join(xc, rc, case.wall_x, case.wall_r)
    case.mode, case.inlet, case.n_ft = "chamber", None, 8.0
    case.breaks = br + [0.0]
    case.ref["rc_up"] = rc_up
    return case


def build_cases() -> list[Case]:
    from rockmoc import (DesignSpec, design_conical, design_ideal, design_mln, design_planar_mln_nasa,
                         design_production, design_rao, design_rao_control, design_rao_toc, design_tic)
    from rockmoc.analysis.boundary_layer import correct_wall_bl
    from rockmoc.thermo import PerfectGas

    G = Gas(1.22, 355.0)
    pg = PerfectGas(1.22, 355.0)

    def spec(Me=3.0, gas=pg, n=40, **k):
        return DesignSpec(gas=gas, throat_radius=0.05, exit_mach=Me, n_characteristics=n, **k)

    cs: list[Case] = []
    add = cs.append
    # ---- Tier 1: inviscid, perfect gas, supersonic inlet at the MOC start line
    add(from_design("t1_mln_M3", "T1", "MLN Me=3 g=1.22", design_mln(spec()), G))
    for Me in (2.0, 4.0, 5.0):
        add(from_design(f"t1_mln_M{Me:g}", "T1", f"MLN Me={Me:g} g=1.22", design_mln(spec(Me)), G))
    for g in (1.15, 1.4):
        gg = PerfectGas(g, 355.0 if g < 1.3 else 287.05)
        add(from_design(f"t1_mln_g{g:g}", "T1", f"MLN Me=3 g={g:g}", design_mln(spec(gas=gg)), Gas(g, gg.R)))
    for n in (10, 20, 80):
        add(from_design(f"t1_mln_n{n}", "T1", f"MLN Me=3, {n} characteristics", design_mln(spec(n=n)), G))
    add(from_design("t1_ideal_M3", "T1", "Ideal (gradual) Me=3", design_ideal(spec()), G))
    for f in (0.6, 0.8):
        add(from_design(f"t1_tic_{int(f*100)}", "T1", f"TIC {int(f*100)}% Me=3", design_tic(spec(truncate_fraction=f)), G))
    for f in (0.6, 0.8, 1.0):
        add(from_design(f"t1_rao_{int(f*100)}", "T1", f"Rao TOP {int(f*100)}% bell Me=3", design_rao(spec(), length_fraction=f), G))
    add(from_design("t1_raotoc_80", "T1", "Rao TOC 80% Me=3", design_rao_toc(spec()), G))
    add(from_design("t1_raoctl_80", "T1", "Rao control-surface (MOC) 80% Me=3", design_rao_control(spec()), G))
    for a in (15.0, 10.0):
        add(from_design(f"t1_cone_{a:g}", "T1", f"Conical {a:g} deg Me=3", design_conical(spec(), half_angle_deg=a), G))
    pr = design_production(spec(n=120), apply_bl=False, kind="mln")
    add(from_design("t1_prod_mln", "T1", "Production MLN (n=120) Me=3", pr.design, G))
    # planar MLN (NASA TM-X-1502 algorithm), h_t = 50 mm, air
    p = design_planar_mln_nasa(2.4, 1.4, store_net=True)
    ht = 0.05
    air = Gas(1.4, 287.05)
    from rockmoc.thermo import PerfectGas as _PG
    pa = _PG(1.4, 287.05)
    nx, ny = np.asarray(p.net_x) * ht, np.asarray(p.net_y) * ht
    nM = np.array([pa.mach_from_prandtl_meyer(v) for v in np.asarray(p.net_nu).ravel()]).reshape(np.shape(p.net_nu))
    ok = np.isfinite(nx) & np.isfinite(nM)
    ref = {"wall_x": (p.contour.x * ht).tolist(), "wall_M": (p.contour.M if p.contour.M is not None else np.full(p.contour.x.size, np.nan)).tolist(),
           "Me_design": 2.4, "kind": "planar_mln", "has_net": False, "uniform_exit": True,
           "Isp_vac": ideal_isp(1.4, 287.05, 300.0, 2.4)}  # stored planar net covers only the corner kernel: not compared
    add(Case("t1_planar_M2.4", "T1", "Planar MLN Me=2.4 g=1.4 (TM-X-1502)", p.contour.x * ht, p.contour.r * ht, ht, air, 1e6, 300.0,
             inlet=isen_inlet(air, 1e6, 300.0, 1.0 + 1e-3 + 0.02), planar=True, ref=ref))
    # ---- GCI (mesh independence) on three representative shapes
    for base in ("t1_mln_M3", "t1_rao_80", "t1_cone_15"):
        b = next(c for c in cs if c.name == base)
        for s, tag in ((0.7071, "coarse"), (1.4142, "fine")):
            c2 = Case(**{**b.__dict__, "name": f"{base}_{tag}", "tier": "GCI", "scale": s, "title": f"{b.title} [{tag}]"})
            cs.append(c2)

    # ---- Tier 2: full nozzle (chamber + converging + throat), transonic / Cd
    from rockmoc.moc.transonic import discharge_coefficient

    for rc in (0.625, 1.0, 2.0):
        c = full_nozzle(from_design(f"t2_rao80_rc{rc:g}", "T2", f"Full Rao 80% nozzle, upstream Rc/Rt={rc:g}",
                                    design_rao(spec(), length_fraction=0.8), G), rc_up=rc)
        c.ref["Cd"] = discharge_coefficient(1.22, rc)
        add(c)
    for nm, res, t in (("t2_ideal_full", design_ideal(spec()), "Full ideal-contour nozzle Me=3 (throat arc Rc/Rt=1.5 both sides)"),
                       ("t2_tic80_full", design_tic(spec(truncate_fraction=0.8)), "Full TIC 80% nozzle Me=3 (throat arc Rc/Rt=1.5)")):
        c = full_nozzle(from_design(nm, "T2", t, res, G), rc_up=1.5)
        c.ref["Cd"] = discharge_coefficient(1.22, 1.5)
        add(c)
    c = full_nozzle(from_design("t2_mln_full", "T2", "Full MLN Me=3 (real converging section, Rc/Rt=1.5)", design_mln(spec()), G))
    c.ref["Cd"] = discharge_coefficient(1.22, 1.5)
    add(c)
    # Cuffel/Back/Massier 45-15 deg, Rc/Rt=0.625, air (geometry from the 1969 paper; chamber CR approximated)
    from rockmoc.validate.experimental import _D5377, _d5377_contour

    Rt_c = 0.8 * 0.0254
    cone = design_conical(DesignSpec(gas=pa, throat_radius=Rt_c, area_ratio=6.63, n_characteristics=20), half_angle_deg=15.0)
    c = full_nozzle(from_design("t2_cuffel", "T2", "Cuffel 1969 45/15 deg cone, Rc/Rt=0.625, air", cone, air, p0=5.2e5, T0=300.0),
                    rc_up=0.625, half_deg=45.0, cr=6.0)
    c.ref.update(Cd=discharge_coefficient(1.4, 0.625), Cd_exp=0.985, M_axis_throat_exp=0.8, M_wall_throat_exp=1.4)
    add(c)
    # NASA TN D-5377 30-15 deg cone with measured wall pressures (inviscid run; viscous below)
    w, Rt5 = _d5377_contour(chamber=True)
    d = _D5377
    ref5 = {"wall_x": w.x.tolist(), "wall_M": w.M.tolist(), "Me_design": float(w.M[-1]), "kind": "d5377", "has_net": False,
            "exp_xd": d["xd"], "exp_pw": d["pw"], "dstar": d["dstar"], "Cd": discharge_coefficient(1.4, d["Rc"] / Rt5)}
    add(Case("t2_d5377_inv", "T2", "NASA TN D-5377 30/15 cone, inviscid", w.x, w.r, Rt5, air, d["P0"], d["T0"], mode="chamber",
             n_ft=8.0, breaks=[float(w.x[np.argmin(np.abs(w.x + d['Rc'] * math.sin(math.radians(30))))]), 0.0], ref=dict(ref5)))

    # ---- Tier 3: real gas (frozen LOX/RP-1 via RocketCEA, JANAF fit in OpenFOAM)
    try:
        cs += tier3(spec)
    except Exception as e:  # noqa: BLE001
        print("tier 3 skipped:", e)

    # ---- Tier 4: viscous (kOmegaSST + Spalding wall function)
    from rockmoc.analysis.boundary_layer import bartz_viscosity
    from rockmoc.validate.experimental import _d5377_bl

    air_v = Gas(1.4, 287.05, sutherland=(1.458e-6, 110.4))
    for inlet in ("adiabatic", "cooled"):
        bl, xd = _d5377_bl(inlet)
        r5 = dict(ref5)
        r5.update(bl_x=bl.corrected.x.tolist() if hasattr(bl, "corrected") else w.x.tolist(), bl_dstar=bl.delta_star.tolist(),
                  bl_theta=bl.theta_mom.tolist(), bl_cf=(bl.cf.tolist() if bl.cf is not None else None),
                  exp_bl={str(k): v for k, v in d["surveys"][inlet].items()})
        Tw = (list(np.array(d["xd"]) * d["dstar"]), d["tw_ad" if inlet == "adiabatic" else "tw_co"])
        add(Case(f"t4_d5377_{inlet}", "T4", f"TN D-5377 viscous, {inlet}-inlet wall temperature", w.x, w.r, Rt5, air_v, d["P0"], d["T0"],
                 mode="chamber", n_ft=1.7, viscous={"Tw": Tw, "dy_wall": 5.0e-6}, nr=100, breaks=ref5_breaks(w, d), ref=r5))  # dy: y+ ~3 downstream, ~12 throat (explicit dt limit)
    M_mol = RU / 355.0
    Gv = Gas(1.22, 355.0, sutherland=sutherland_fit(lambda T: bartz_viscosity(T, M_mol)))
    prv = design_production(spec(n=120), apply_bl=True, kind="mln", pc=1e6, Tc=3200.0)
    c = from_design("t4_prod_bl", "T4", "Production MLN with BL correction, viscous, pc=1 MPa", prv.design, Gv, p0=1e6)
    c.wall_x, c.wall_r = prv.bl.corrected.x.copy(), prv.bl.corrected.r.copy()
    perf = prv.field.summary()  # v0.10: exit-plane BL control volume, not wall-shear integral
    c.ref.update(Isp_vac=perf["Isp_vac_s"], Cf_vac=perf["Cf_vac"], bl_x=prv.bl.inviscid.x.tolist(), bl_dstar=prv.bl.delta_star.tolist(),
                 bl_theta=prv.bl.theta_mom.tolist(), bl_cf=prv.bl.cf.tolist() if prv.bl.cf is not None else None,
                 dcf_friction=prv.bl.meta["delta_cf_friction"])
    c.viscous, c.nr, c.n_ft = {"Tw": None, "dy_wall": 5.0e-6}, 100, 2.0  # y+ ~1 at throat (pc 1 MPa)
    add(c)
    c_mln = design_mln(spec(n=120))
    c = from_design("t4_mln_inviscid_wall_visc", "T4", "MLN (no BL correction) viscous, pc=1 MPa: cost of skipping BL", c_mln, Gv, p0=1e6)
    c.viscous, c.nr, c.n_ft = {"Tw": None, "dy_wall": 5.0e-6}, 100, 2.0  # y+ ~1 at throat (pc 1 MPa)
    blu = correct_wall_bl(c_mln.contour, pg, throat_radius=0.05, pc=1e6, Tc=3200.0, turbulent=True)
    perf = c_mln.performance(pc=1e6, Tc=3200.0, delta_cf_friction=blu.meta["delta_cf_uncorrected_vac"])
    c.ref.update(Isp_vac=perf["Isp_vac_s"], Cf_vac=perf["Cf_vac"])
    add(c)
    # sea-level overexpanded: flow separation vs Schmucker
    from rockmoc.analysis.separation import separated_cf

    for nm, res, pc, conical in (("rao80_pc1", design_rao(spec(), length_fraction=0.8), 1.0e6, False),
                                 ("rao80_pc1.5", design_rao(spec(), length_fraction=0.8), 1.5e6, False),
                                 ("cone15_pc1", design_conical(spec(), half_angle_deg=15.0), 1.0e6, True)):
        c = from_design(f"t4_sep_{nm}", "T4", f"Sea-level separation {nm}", res, Gv, p0=pc)
        # rockmoc v0.10 path: MOC wall pressure locates separation, Kalt-Badal recovery, losses on vacuum part
        pf = res.performance(pc=pc, Tc=3200.0, pa=101325.0, separation="schmucker")
        sm = pf["separation"]
        att = res.performance(pc=pc, Tc=3200.0, pa=101325.0)
        c.ref.update(sep_separated=sm["separated"], sep_eps=sm["eps_sep"], sep_Cf=pf["Cf_amb"], sep_Cf_att=att["Cf_amb"],
                     sep_psep=sm["psep_pa"], sep_wall_pressure=sm["wall_pressure"])
        sr = separated_cf(pg, area_ratio=float((res.contour.r[-1] / 0.05) ** 2), pc=pc, pa=101325.0, conical=conical)
        c.ref.update(sep_eps_v09=sr.eps_sep, sep_Cf_v09=sr.Cf)  # old 1-D model, for the report
        c.viscous, c.nr, c.n_ft = {"Tw": None, "dy_wall": 1.0e-5}, 100, 3.0  # y+ ~1-2 downstream, ~4 throat
        c.plume = {"pa": 101325.0, "Ta": 300.0}
        add(c)
    return cs


def ref5_breaks(w, d):
    return [float(-d["Rc"] * math.sin(math.radians(30))), 0.0]


def tier3(spec):
    from rockmoc import design_mln
    from rockmoc.thermo import PerfectGas, VariableGammaGas
    from rockmoc.thermo.cea_table import build_cea_table

    pc, of = 5e6, 2.6
    out = []
    tbs = {m: build_cea_table("LOX", "RP1", of, pc, mode=m, eps_max=60.0) for m in ("frozen", "equilibrium")}
    tb = tbs["frozen"]
    Tc = tb.Tc
    rhoc = pc / (tb.R * Tc)
    M, a, T = tb.M, tb.a, tb.T_Tc * Tc
    p = tb.p_pc * pc
    rho = tb.rho_rhoc * rhoc
    Rf = float(np.median(p / (rho * T)))
    V = M * a
    h = -0.5 * V**2
    cp = np.gradient(h, T)
    co = np.polyfit(T, cp / Rf, 4)[::-1]  # a0..a4
    jan = [float(c) for c in co] + [0.0, 0.0]
    g_t = float(np.interp(1.05, tb.M, tb.gamma_local))
    table = {"M": M.tolist(), "p_pc": tb.p_pc.tolist(), "T_Tc": tb.T_Tc.tolist(), "rho_rhoc": tb.rho_rhoc.tolist(), "a": a.tolist(),
             "Tc": Tc, "gamma_c": float(tb.gamma_local[0]), "rho_c": rhoc}
    Gf = Gas(g_t, Rf, kind="janaf", janaf=jan, table=table)
    fit_err = float(np.max(np.abs(np.polyval(co[::-1], T) * Rf / cp - 1)))
    for nm, tbl, title in (("frozen", tb, "MLN designed with frozen CEA table, frozen-gas CFD"),
                           ("equil", tbs["equilibrium"], "MLN designed with equilibrium CEA table, frozen-gas CFD (bound)")):
        res = design_mln(spec(gas=VariableGammaGas(tbl)))
        c = from_design(f"t3_{nm}", "T3", title, res, Gf, p0=pc, T0=Tc)
        c.ref.update(janaf_fit_max_err=fit_err, gas_table=nm)
        out.append(c)
    gp = PerfectGas(g_t, Rf)
    c = from_design("t3_perfect_contour", "T3", f"MLN designed with perfect gas g={g_t:.3f}, frozen-gas CFD", design_mln(spec(gas=gp)), Gf, p0=pc, T0=Tc)
    c.ref["gas_table"] = "perfect"
    out.append(c)
    return out


# ----------------------------------------------------------------------------- mesh
class Block:
    def __init__(s, planar):
        s.planar, s.V, s.key, s.blocks, s.edges, s.faces = planar, [], {}, [], [], {}
        h = math.radians(WEDGE / 2)
        s.ch, s.sh = math.cos(h), math.sin(h)

    def v(s, x, r, side):
        if not s.planar and abs(r) < 1e-14:
            k = (round(x, 12), 0.0, "a")
            p = (x, 0.0, 0.0)
        else:
            k = (round(x, 12), round(r, 12), side)
            p = (x, r, (-1 if side == 0 else 1) * 0.5 * s.dz) if s.planar else (x, r * s.ch, (-1 if side == 0 else 1) * r * s.sh)
        if k not in s.key:
            s.key[k] = len(s.V)
            s.V.append(p)
        return s.key[k]

    def pt(s, x, r, side):
        if s.planar:
            return (x, r, (-1 if side == 0 else 1) * 0.5 * s.dz)
        return (x, r * s.ch, (-1 if side == 0 else 1) * r * s.sh)

    def quad(s, x0, x1, lo0, lo1, hi0, hi1, n, grading, patches, top_pts=None):
        hx = [s.v(x0, lo0, 0), s.v(x1, lo1, 0), s.v(x1, hi1, 0), s.v(x0, hi0, 0),
              s.v(x0, lo0, 1), s.v(x1, lo1, 1), s.v(x1, hi1, 1), s.v(x0, hi0, 1)]
        gr = grading[1] if isinstance(grading[1], str) else f"{grading[1]:.6g}"
        s.blocks.append(f"hex ({' '.join(map(str, hx))}) ({n[0]} {n[1]} 1) simpleGrading ({grading[0]:.6g} {gr} 1)")
        if top_pts is not None and len(top_pts):
            for side, (a, b) in ((0, (3, 2)), (1, (7, 6))):
                pts = " ".join(f"({p[0]:.9g} {p[1]:.9g} {p[2]:.9g})" for p in (s.pt(x, r, side) for x, r in top_pts))
                s.edges.append(f"polyLine {hx[a]} {hx[b]} ({pts})")
        fdef = {"xmin": (0, 4, 7, 3), "xmax": (1, 2, 6, 5), "rmin": (0, 1, 5, 4), "rmax": (3, 7, 6, 2),
                "back": (0, 3, 2, 1), "front": (4, 5, 6, 7)}
        for side, patch in patches.items():
            s.faces.setdefault(patch, []).append("(" + " ".join(str(hx[i]) for i in fdef[side]) + ")")

    def write(s, path, ptypes):
        V = "\n".join(f"    ({x:.9g} {y:.9g} {z:.9g})" for x, y, z in s.V)
        B = "\n".join("    " + b for b in s.blocks)
        E = "\n".join("    " + e for e in s.edges)
        P = "\n".join(f"    {n} {{ type {ptypes[n]}; faces ({' '.join(f)}); }}" for n, f in s.faces.items())
        wr(path, f"{hdr('dictionary', 'blockMeshDict')}scale 1;\nvertices\n(\n{V}\n);\nblocks\n(\n{B}\n);\n"
                 f"edges\n(\n{E}\n);\nboundary\n(\n{P}\n);\nmergePatchPairs ();\n")


def grading_for(dy_frac, n):
    """Expansion ratio last/first (axis->wall) whose last cell is dy_frac of the height."""
    lo, hi = 0.5, 1.0
    for _ in range(100):
        q = 0.5 * (lo + hi)
        last = (q - 1) * q ** (n - 1) / (q**n - 1)
        lo, hi = (q, hi) if last < dy_frac else (lo, q)
    return q ** (n - 1)


def smooth_wall(x, r, Rt, n=500, n_throat=200):
    """PCHIP resample in arc length (C1, passes through every given point, no overshoot at
    sharp cone/MLN corners), densified around the throat. A raw polyline of a coarse
    contour puts slope kinks into the mesh wall, each of which launches a wave."""
    from scipy.interpolate import PchipInterpolator

    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(r)))])
    st = float(np.interp(0.0, x, s)) if x[0] < 0 < x[-1] else 0.0
    sn = np.unique(np.concatenate([s, np.linspace(0, s[-1], n),
                                   np.linspace(max(st - 1.5 * Rt, 0), min(st + 2.0 * Rt, s[-1]), n_throat)]))
    return PchipInterpolator(s, x)(sn), PchipInterpolator(s, r)(sn)


def radial_grading(case, nr):
    """Inviscid: mild 2:1 axis->wall. Viscous: uniform core (75 % height, half the cells) +
    geometric wall zone whose last cell is dy_wall at the throat (relative grading, so all
    axial blocks share it; dy grows with r downstream where y+ is lower anyway)."""
    if not case.viscous:
        return 0.5
    nw = nr // 2
    gw = grading_for(case.viscous["dy_wall"] / (0.25 * case.Rt), nw)
    return f"((0.75 {nr - nw} 1) (0.25 {nw} {gw:.6g}))"


def mesh(case: Case, path: Path):
    b = Block(case.planar)
    b.dz = 0.1 * case.Rt
    x, r = smooth_wall(case.wall_x, case.wall_r, case.Rt)
    nr = max(8, round(case.nr * case.scale))
    g = radial_grading(case, nr)
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(r)))])
    K = 30
    idx = {0, len(x) - 1} | {int(np.argmin(np.abs(s - t))) for t in np.linspace(0, s[-1], K + 1)}
    idx |= {int(np.argmin(np.abs(x - xb))) for xb in case.breaks}
    idx |= {int(np.argmin(np.abs(x - xb * case.Rt))) for xb in (-1.0, -0.5, 0.5, 1.0, 2.0) if x[0] < xb * case.Rt < x[-1]}
    idx = sorted(idx)
    keep = [idx[0]]
    for i in idx[1:]:
        if i - keep[-1] >= 2 or i == len(x) - 1:
            keep.append(i)
    if keep[-1] - keep[-2] < 2 and len(keep) > 2:
        keep.pop(-2)
    nx_tot = 0
    for k, (i0, i1) in enumerate(zip(keep[:-1], keep[1:])):
        L = x[i1] - x[i0]
        rm = 0.5 * (r[i0] + r[i1])
        xm = 0.5 * (x[i0] + x[i1]) / case.Rt
        nx = max(2, round(L / (case.ar * rm / nr) * (1.0 + 2.0 * math.exp(-xm * xm))))  # 3x axial density at the throat
        nx_tot += nx
        pat = {"rmax": "wall", "rmin": "axis", "back": "back", "front": "front"}
        if k == 0:
            pat["xmin"] = "inlet"
        if i1 == len(x) - 1 and not case.plume:
            pat["xmax"] = "outlet"
        b.quad(x[i0], x[i1], 0.0, 0.0, r[i0], r[i1], (nx, nr), (1.0, g), pat, list(zip(x[i0 + 1:i1], r[i0 + 1:i1])))
    ncell = nx_tot * nr
    if case.plume:
        xe, Re = x[-1], r[-1]
        Lp, Rf = 5.0 * Re, 3.5 * Re
        dx_e = case.ar * Re / nr
        nxp = max(20, round(case.scale * 70))
        gx = grading_ratio(Lp, nxp, dx_e)
        no = max(10, round(case.scale * 40))
        if case.viscous:
            # free shear layer, not a wall: starting the ambient block at the wall-cell size
            # (2.5e-5 m) made viscous+plume runs diverge in the first steps; core-cell size is stable
            d_last = Re / nr
        else:
            q = g ** (1 / (nr - 1))
            d_last = Re * (1 - q) / (1 - q**nr) * g
        go = grading_ratio(Rf - Re, no, d_last)
        b.quad(xe, xe + Lp, 0.0, 0.0, Re, Re, (nxp, nr), (gx, g), {"xmax": "outlet", "rmin": "axis", "back": "back", "front": "front"})
        b.quad(xe, xe + Lp, Re, Re, Rf, Rf, (nxp, no), (gx, go),
               {"xmin": "ambient", "xmax": "outlet", "rmax": "ambient", "back": "back", "front": "front"})
        ncell += nxp * (nr + no)
    ptypes = {"inlet": "patch", "outlet": "patch", "wall": "wall", "ambient": "patch",
              "axis": "symmetryPlane" if case.planar else "empty",
              "front": "empty" if case.planar else "wedge", "back": "empty" if case.planar else "wedge"}
    b.write(path / "system" / "blockMeshDict", ptypes)
    return ncell, nr, g


def grading_ratio(L, n, d0):
    """Expansion ratio last/first for n cells over L with first cell d0."""
    if d0 * n >= L:
        return 1.0
    lo, hi = 1.0 + 1e-9, 2.0
    for _ in range(100):
        q = 0.5 * (lo + hi)
        if d0 * (q**n - 1) / (q - 1) < L:
            lo = q
        else:
            hi = q
    return q ** (n - 1)


# ----------------------------------------------------------------------------- OpenFOAM case
def hdr(cls, obj):
    return f"FoamFile\n{{\n    version 2.0; format ascii; class {cls}; object {obj};\n}}\n"


def wr(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.open("w", newline="\n").write(text)


def bc_block(case: Case, entries: dict):
    common = {"axis": "type symmetryPlane;" if case.planar else "type empty;",
              "front": "type empty;" if case.planar else "type wedge;",
              "back": "type empty;" if case.planar else "type wedge;"}
    pats = ["inlet", "outlet", "wall"] + (["ambient"] if case.plume else [])
    body = "\n".join(f"    {p} {{ {entries[p]} }}" for p in pats)
    body += "\n" + "\n".join(f"    {p} {{ {v} }}" for p, v in common.items())
    return f"boundaryField\n{{\n{body}\n}}\n"


def field_file(case, cls, name, dim, init, entries):
    return hdr(cls, name) + f"dimensions {dim};\ninternalField uniform {init};\n" + bc_block(case, entries)


def write_case(case: Case) -> dict:
    d = CASES_DIR / case.name
    if d.exists():
        shutil.rmtree(d)
    (d / "system").mkdir(parents=True)
    ncell, nr, g = mesh(case, d)
    gas = case.gas
    W = RU / gas.R
    pa = case.plume["pa"] if case.plume else None
    Ta = case.plume["Ta"] if case.plume else 300.0
    vis = case.viscous is not None
    # thermo
    if gas.kind == "janaf":
        c = " ".join(f"{v:.10g}" for v in gas.janaf)
        thermo, td = "janaf", f"Tlow 200; Thigh 6000; Tcommon 1000; highCpCoeffs ({c}); lowCpCoeffs ({c});"
    else:
        thermo, td = "hConst", f"Cp {gas.gamma * gas.R / (gas.gamma - 1):.10g}; Hf 0;"
    if vis or gas.kind == "janaf":  # janaf is only compiled with sutherland; As=0 -> inviscid
        As, Ts = gas.sutherland if vis else (0.0, 1.0)
        transport, tr = "sutherland", f"As {As:.6g}; Ts {Ts:.6g};"
    else:
        transport, tr = "const", "mu 0; Pr 1;"
    wr(d / "constant/thermophysicalProperties", hdr("dictionary", "thermophysicalProperties") +
       f"thermoType {{ type hePsiThermo; mixture pureMixture; transport {transport}; thermo {thermo}; equationOfState perfectGas; "
       f"specie specie; energy sensibleInternalEnergy; }}\nmixture {{ specie {{ molWeight {W:.10g}; }} thermodynamics {{ {td} }} "
       f"transport {{ {tr} }} }}\n")
    turb = "simulationType RAS;\nRAS { model kOmegaSST; RASModel kOmegaSST; turbulence on; printCoeffs on; }\n" if vis else "simulationType laminar;\n"
    wr(d / "constant/turbulenceProperties", hdr("dictionary", "turbulenceProperties") + turb)
    wr(d / "constant/momentumTransport", hdr("dictionary", "momentumTransport") + turb)
    # 1-D speeds for run length
    g_ = gas.gamma
    xs = np.linspace(case.wall_x[0], case.wall_x[-1], 400)
    rs = np.interp(xs, case.wall_x, case.wall_r)
    T0 = case.T0 if case.mode == "chamber" else case.inlet[2] * (1 + 0.5 * (g_ - 1) * case.inlet[0] ** 2)
    a0 = math.sqrt(g_ * gas.R * T0)
    # rough 1-D speed: subsonic part ~ continuity from the throat, supersonic ~0.9 a0
    # (chamber speed floored at 0.3 a0: the 1-D initial field already holds the chamber state,
    # only acoustic settling is needed there)
    u = np.where(xs < 0, np.maximum(0.6 * a0 * (case.Rt / rs) ** 2, 0.3 * a0), 0.9 * a0)
    t_ft = float(np.sum(np.diff(xs) / u[:-1]))
    end = case.n_ft * t_ft
    if case.plume:
        end += 4.0 * 5.0 * case.wall_r[-1] / (0.9 * a0)
    # fields
    if case.mode == "throat":
        Mi, pi, Ti = case.inlet
        Ui = Mi * math.sqrt(g_ * gas.R * Ti)
        inl = {"p": f"type fixedValue; value uniform {pi:.8g};", "T": f"type fixedValue; value uniform {Ti:.8g};",
               "U": f"type fixedValue; value uniform ({Ui:.8g} 0 0);"}
    else:
        inl = {"p": f"type totalPressure; p0 uniform {case.p0:.8g}; value uniform {case.p0:.8g};",
               "T": f"type totalTemperature; gamma {g_}; T0 uniform {case.T0:.8g}; value uniform {case.T0:.8g};",
               "U": "type pressureInletOutletVelocity; value uniform (0 0 0);"}
    if case.plume:
        outl = {"p": f"type waveTransmissive; field p; phi phi; rho rho; psi thermo:psi; gamma {g_}; fieldInf {pa}; lInf {5 * case.wall_r[-1]:.6g}; value uniform {pa};",
                "T": f"type inletOutlet; inletValue uniform {Ta}; value uniform {Ta};", "U": "type zeroGradient;"}
        amb = {"p": f"type totalPressure; p0 uniform {pa}; value uniform {pa};", "T": f"type inletOutlet; inletValue uniform {Ta}; value uniform {Ta};",
               "U": "type pressureInletOutletVelocity; value uniform (0 0 0);"}
    else:
        outl = {k: "type zeroGradient;" for k in ("p", "T", "U")}
        amb = {}
    Tw = case.viscous.get("Tw") if vis else None
    if Tw is None:
        wT = "type zeroGradient;"
    elif isinstance(Tw, (int, float)):
        wT = f"type fixedValue; value uniform {Tw};"
    else:
        wT = "type fixedValue; value uniform TWALL;"
    wall = {"p": "type zeroGradient;", "T": wT, "U": "type noSlip;" if vis else "type slip;"}

    def E(f):
        e = {"inlet": inl[f], "outlet": outl[f], "wall": wall[f]}
        if case.plume:
            e["ambient"] = amb[f]
        return e

    wr(d / "0/p", field_file(case, "volScalarField", "p", "[1 -1 -2 0 0 0 0]", 1e5, E("p")))
    wr(d / "0/T", field_file(case, "volScalarField", "T", "[0 0 0 1 0 0 0]", 300, E("T")))
    wr(d / "0/U", field_file(case, "volVectorField", "U", "[0 1 -1 0 0 0 0]", "(0 0 0)", E("U")))
    if vis:
        k0, om0 = 100.0, 1e4
        L = 0.07 * case.wall_r[0]
        tk = {"inlet": "type turbulentIntensityKineticEnergyInlet; intensity 0.05; value uniform 100;",
              "outlet": "type inletOutlet; inletValue uniform 1; value uniform 1;", "wall": "type kqRWallFunction; value uniform 100;",
              "ambient": "type inletOutlet; inletValue uniform 1; value uniform 1;"}
        to = {"inlet": f"type turbulentMixingLengthFrequencyInlet; mixingLength {L:.6g}; value uniform 1e4;",
              "outlet": "type inletOutlet; inletValue uniform 100; value uniform 100;", "wall": "type omegaWallFunction; value uniform 1e4;",
              "ambient": "type inletOutlet; inletValue uniform 100; value uniform 100;"}
        tn = {"inlet": "type calculated; value uniform 0;", "outlet": "type calculated; value uniform 0;",
              "wall": "type nutUSpaldingWallFunction; value uniform 0;", "ambient": "type calculated; value uniform 0;"}
        ta = {"inlet": "type calculated; value uniform 0;", "outlet": "type calculated; value uniform 0;",
              "wall": "type compressible::alphatWallFunction; Prt 0.85; value uniform 0;", "ambient": "type calculated; value uniform 0;"}
        wr(d / "0/k", field_file(case, "volScalarField", "k", "[0 2 -2 0 0 0 0]", k0, tk))
        wr(d / "0/omega", field_file(case, "volScalarField", "omega", "[0 0 -1 0 0 0 0]", om0, to))
        wr(d / "0/nut", field_file(case, "volScalarField", "nut", "[0 2 -1 0 0 0 0]", 0, tn))
        wr(d / "0/alphat", field_file(case, "volScalarField", "alphat", "[1 -1 -1 0 0 0 0]", 0, ta))
    # system
    fo = """    flowRates { type surfaceFieldValue; libs (fieldFunctionObjects); writeControl timeStep; writeInterval 50; log false; writeFields false; regionType patch; name inlet; operation sum; fields (phi); }
    flowRateOut { $flowRates; name outlet; }
    machNo { type MachNo; libs (fieldFunctionObjects); writeControl writeTime; }
"""
    if case.mode == "chamber":  # chamber acoustics (central scheme at low M): time-average the last 40 %
        fo += (f"    fieldAverage1 {{ type fieldAverage; libs (fieldFunctionObjects); writeControl writeTime; "
               f"timeStart {0.6 * end:.6g}; fields ( U {{ mean on; prime2Mean off; base time; }} "
               f"p {{ mean on; prime2Mean off; base time; }} T {{ mean on; prime2Mean off; base time; }} "
               f"rho {{ mean on; prime2Mean off; base time; }} ); }}\n")
    if vis:
        fo += """    wallShearStress { type wallShearStress; libs (fieldFunctionObjects); patches (wall); writeControl writeTime; log false; }
    yPlus { type yPlus; libs (fieldFunctionObjects); writeControl writeTime; log false; }
"""
    wr(d / "system/controlDict", hdr("dictionary", "controlDict") +
       f"application rhoCentralFoamAxi;\nstartFrom latestTime;\nstartTime 0;\nstopAt endTime;\nendTime {end:.6g};\ndeltaT 1e-10;\n"
       f"writeControl adjustable;\nwriteInterval {end / 4:.6g};\npurgeWrite 2;\nwriteFormat ascii;\nwritePrecision 8;\n"
       f"timeFormat general;\ntimePrecision 8;\nrunTimeModifiable true;\nadjustTimeStep yes;\nmaxCo 0.4;\nmaxDeltaT 1;\nTmin 20;\n"
       f"functions\n{{\n{fo}}}\n")
    turb_div = "div(phi,k) Gauss limitedLinear 1; div(phi,omega) Gauss limitedLinear 1;" if vis else ""
    wr(d / "system/fvSchemes", hdr("dictionary", "fvSchemes") +
       "fluxScheme Kurganov;\nddtSchemes { default Euler; }\ngradSchemes { default Gauss linear; }\n"
       f"divSchemes {{ default none; div(tauMC) Gauss linear; {turb_div} }}\nlaplacianSchemes {{ default Gauss linear corrected; }}\n"
       "interpolationSchemes { default linear; reconstruct(rho) vanLeer; reconstruct(U) vanLeerV; reconstruct(T) vanLeer; }\n"
       "snGradSchemes { default corrected; }\nwallDist { method meshWave; }\n")
    wr(d / "system/fvSolution", hdr("dictionary", "fvSolution") +
       'solvers\n{\n    "(rho|rhoU|rhoE)" { solver diagonal; }\n    U { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-09; relTol 0.01; }\n'
       '    e { $U; tolerance 1e-10; relTol 0; }\n    "(k|omega)" { $U; tolerance 1e-10; relTol 0; }\n}\n')
    wr(d / "system/decomposeParDict", hdr("dictionary", "decomposeParDict") + f"numberOfSubdomains {NP};\nmethod scotch;\n")
    meta = {"name": case.name, "planar": case.planar, "mode": case.mode, "wall_x": case.wall_x.tolist(), "wall_r": case.wall_r.tolist(),
            "Rt": case.Rt, "gamma_init": g_, "R_init": gas.R, "p0": case.p0, "T0": case.T0, "plume": case.plume,
            "inlet_M": case.inlet[0] if case.inlet else None, "inlet_p": case.inlet[1] if case.inlet else None,
            "inlet_T": case.inlet[2] if case.inlet else None, "Tw_profile": Tw if (vis and isinstance(Tw, (list, tuple))) else None,
            "ncell": ncell, "endTime": end}
    wr(d / "meta.json", json.dumps(meta))
    ref = {**case.ref, "title": case.title, "tier": case.tier, "Rt": case.Rt, "p0": case.p0, "T0": case.T0, "gamma": gas.gamma,
           "R": gas.R, "table": gas.table, "planar": case.planar, "mode": case.mode, "viscous": vis, "plume": case.plume,
           "ncell": ncell, "nr": nr, "scale": case.scale, "wall_geom_x": case.wall_x.tolist(), "wall_geom_r": case.wall_r.tolist()}
    wr(d / "ref.json", json.dumps(ref, default=float))
    wr(d / "Allrun", ALLRUN)
    return meta


ALLRUN = r"""#!/bin/bash
# usage: bash Allrun <cpu-set>
source /usr/lib/openfoam/openfoam2412/etc/bashrc || true
cd "$(dirname "$0")"
[ -f DONE ] && exit 0
rm -rf processor* [1-9]* 0.[0-9]* FAILED
fail() { echo "$1" > FAILED; echo "$(date +%T) $(basename $PWD) FAILED $1" >> ../queue.log; exit 1; }
echo "$(date +%T) $(basename $PWD) start" >> ../queue.log
blockMesh > log.blockMesh 2>&1 || fail blockMesh
checkMesh > log.checkMesh 2>&1
postProcess -func writeCellCentres -time 0 > log.cc 2>&1 || fail cellCentres
~/venv/bin/python ../init_fields.py . > log.init 2>&1 || fail init
rm -f 0/C 0/Cx 0/Cy 0/Cz 0/V
decomposePar -force > log.decompose 2>&1 || fail decompose
mpirun --bind-to none -np 8 rhoCentralFoamAxi -parallel > log.solver 2>&1 || fail solver
reconstructPar -latestTime > log.reconstruct 2>&1 || fail reconstruct
~/venv/bin/python ../extract.py . > log.extract 2>&1 || fail extract
rm -rf processor*
touch DONE
echo "$(date +%T) $(basename $PWD) done" >> ../queue.log
"""


# ----------------------------------------------------------------------------- remote
def lane_script(cases, cpus, lane):
    lst = " ".join(cases)
    return f"""#!/bin/bash
cd {REMOTE}
for c in {lst}; do bash $c/Allrun {cpus}; done
echo "$(date +%T) lane {lane} finished" >> queue.log
"""


def launch(names):
    import tarfile
    import tempfile

    import remote_run as rr

    order = names
    # greedy balance by cells * steps proxy (stored in meta)
    cost = {}
    for n in order:
        m = json.loads((CASES_DIR / n / "meta.json").read_text())
        cost[n] = m["ncell"] * (3.0 if json.loads((CASES_DIR / n / "ref.json").read_text())["viscous"] else 1.0) * (m["ncell"] ** 0.5)
    lanes = {"A": [], "B": []}
    load = {"A": 0.0, "B": 0.0}
    for n in order:  # keep given priority order, assign to lighter lane
        k = min(load, key=load.get)
        lanes[k].append(n)
        load[k] += cost[n]
    tgz = Path(tempfile.gettempdir()) / "cfdval.tgz"
    with tarfile.open(tgz, "w:gz") as t:
        for n in order:
            t.add(CASES_DIR / n, arcname=f"cfdval/{n}")
        for f in ("init_fields.py", "extract.py"):
            t.add(HERE / "remote" / f, arcname=f"cfdval/{f}")
    rr.scp(str(tgz), f"{rr.HOST}:{rr.XFER_WIN}/cfdval.tgz")
    rr.wsl(f"mkdir -p ~/foamstudio/nozzles && cd ~/foamstudio/nozzles && tar -xzf {rr.XFER_WSL}/cfdval.tgz && echo ok")
    for k, cs in lanes.items():
        if not cs:
            continue
        p = Path(tempfile.gettempdir()) / f"lane{k}.sh"
        p.open("w", newline="\n").write(lane_script(cs, LANES[k], k))
        rr.wsl(f"mkdir -p {REMOTE}/lane{k}")
        rr.run(f"cfdval/lane{k}", str(p))
        print(f"lane {k}: {' '.join(cs)}")


def status():
    import remote_run as rr

    rr.scp(str(HERE / "remote" / "status.sh"), f"{rr.HOST}:{rr.XFER_WIN}/status.sh")
    print(rr.wsl(f"cp {rr.XFER_WSL}/status.sh {REMOTE}/ && bash {REMOTE}/status.sh"))


def pull():
    import tarfile

    import remote_run as rr

    rr.wsl(f"cd {REMOTE} && tar -czf {rr.XFER_WSL}/cfdval_out.tgz $(ls -d */DONE */FAILED 2>/dev/null | xargs -n1 dirname | sort -u | "
           "while read d; do ls $d/cfd.npz $d/log.* $d/FAILED 2>/dev/null; done) queue.log")
    dst = RES_DIR / "raw"
    dst.mkdir(parents=True, exist_ok=True)
    rr.scp(f"{rr.HOST}:{rr.XFER_WIN}/cfdval_out.tgz", str(dst / "out.tgz"))
    with tarfile.open(dst / "out.tgz") as t:
        t.extractall(dst, filter="data")
    print("pulled:", sorted(p.parent.name for p in dst.glob("*/cfd.npz")))


# ----------------------------------------------------------------------------- main
def main():
    a = sys.argv[1:]
    cmd, sel = (a[0] if a else ""), a[1:]
    if cmd == "build":
        cs = build_cases()
        for c in cs:
            if sel and c.name not in sel:
                continue
            m = write_case(c)
            print(f"{c.name:28s} {c.tier:4s} cells={m['ncell']:7d} end={m['endTime']:.3g}s  {c.title}")
    elif cmd == "launch":
        names = sel or [p.name for p in sorted(CASES_DIR.iterdir()) if (p / "meta.json").exists()]
        launch(names)
    elif cmd == "status":
        status()
    elif cmd == "pull":
        pull()
    elif cmd == "compare":
        from compare import compare_all

        compare_all(sel)
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
