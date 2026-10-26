"""Validation against published experiments and independent reference designs.

Every number below is transcribed from the cited document (local copies in
``refs/``). Run :func:`run_experimental_validation` or
``python -m rockmoc.validate.experimental`` to regenerate
``docs/validation_experimental.md``.

Sources
-------
TN-D-5377  Boldman, Schmidt & Ehlers, NASA TN D-5377 (1969): turbulent BL
           surveys in a 30°–15° water-cooled conical nozzle, air, M 2.1–4.4.
CUFFEL69   Cuffel, Back & Massier, AIAA J. 7(7) (1969): transonic flow,
           45°–15° cone, Rc/Rt = 0.625.
TM-107285  Arrington, Reed & Rivera, NASA TM-107285 (1996): wall points of a
           Rao-code (fixed envelope) bell, Table I.
SP-8120    NASA SP-8120 (1976) Rao θ_N/θ_e chart, 80 % bell, as reproduced in
           Newlands, Aspire Space TOP note (2017); read to ±1°.
ANGELINO   Angelino (1964) centred-fan plug contour, derived here from mass
           conservation across each Mach line from the lip.
ZEBBICHE   T. Zebbiche, Doctorat d'État thesis, Univ. Blida (open PDF),
           axisymmetric MLN perfect-gas tables.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field

import numpy as np

from rockmoc.exceptions import MOCFallbackWarning
from rockmoc.geometry.contour import Contour
from rockmoc.thermo.perfect_gas import PerfectGas

IN = 0.0254
AIR = PerfectGas(1.4, 287.05)


@dataclass
class Check:
    nozzle: str
    case: str
    quantity: str
    reference: float
    predicted: float
    tol: float  # absolute (``rel=False``) or relative
    source: str
    rel: bool = True
    note: str = ""
    informational: bool = False
    passed: bool = field(init=False)

    def __post_init__(self) -> None:
        err = abs(self.predicted - self.reference)
        self.passed = err <= (self.tol * abs(self.reference) if self.rel else self.tol)

    @property
    def error(self) -> str:
        if self.rel:
            return f"{100.0 * (self.predicted / self.reference - 1.0):+.1f} %"
        return f"{self.predicted - self.reference:+.3g}"


def _sutherland_air(T):
    T = np.asarray(T, dtype=np.float64)
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


# --- NASA TN D-5377 ------------------------------------------------------------
_D5377 = {
    "dstar": 1.492 * IN, "Rc": 1.5 * IN, "P0": 300 * 6894.757, "T0": 539.0,
    # Table I: x/d*, Pw/P0, Tw adiabatic inlet [K], Tw cooled inlet [K] (inlet z = 0 → throat → exit)
    "xd": [-3.175, -3.026, -2.356, -1.684, -1.446, -1.214, -.979, -.744, -.411, -.117,
           0.0, .087, .171, .263, .425, .818, 1.834, 3.665, 5.497, 6.381, 7.394],
    "pw": [.99932, .99946, .99847, .99565, .99325, .98887, .9804, .9624, .8399, .5733,
           .4600, .3733, .2785, .2084, .1820, .1125, .0375, .0094, .0038, None, .0024],
    "tw_ad": [389, 390, 396, 424, 433, 441, 448, 458, 471, 466, 462, 458, 455, 449, 444, 431, 403, 345, 321, 317, 316],
    "tw_co": [374, 376, 389, 420, 424, 428, 435, 444, 457, 457, 458, 446, 441, 437, 435, 425, 398, 351, 322, 318, 317],
    # Table I(b): nozzle-entrance (θ, δ, Δ) [m]
    "inlet": {"adiabatic": (0.0871e-2, 0.848e-2, 0.00848e-2), "cooled": (0.160e-2, 1.549e-2, 1.402e-2)},
    # Tables II–VII: station x/d*, u_inf [m/s], (θ, δ*) [cm] integrated to edge and to free stream
    "surveys": {
        "adiabatic": {0.818: (706, .0036, .0076, .0086, .0102), 3.665: (891, .013, .014, .061, .066),
                      6.381: (927, .0229, .025, .142, .155)},
        "cooled": {0.818: (708, .0061, .0099, .0048, .0041), 3.665: (888, .0185, .028, .056, .053),
                   6.381: (933, .028, .041, .155, .155)},
    },
}


def _d5377_contour(chamber: bool = True) -> tuple[Contour, float]:
    """30°–15° nozzle wall from the inlet (z = 0, d = 6.5 in) or from the throat."""
    d = _D5377
    Rt, Rc = 0.5 * d["dstar"], d["Rc"]
    a, c = math.radians(15.0), math.radians(30.0)
    x_t = Rc * math.sin(a)
    xc, rc = -Rc * math.sin(c), Rt + Rc * (1 - math.cos(c))
    x0 = xc - (3.25 * IN - rc) / math.tan(c) if chamber else 0.0
    x = np.linspace(x0, 7.394 * d["dstar"], 1500 if chamber else 800)
    r = np.where(x < xc, rc + (xc - x) * math.tan(c),
                 np.where(x <= x_t, Rt + Rc - np.sqrt(np.maximum(Rc**2 - x**2, 0.0)),
                          Rt + Rc * (1 - math.cos(a)) + (x - x_t) * math.tan(a)))
    th = np.where(x < xc, -c, np.where(x <= x_t, np.arcsin(np.clip(x / Rc, -1.0, 1.0)), a))
    # Edge pressure: Table I wall taps + survey free-stream Mach (u_inf, T0). The
    # throat tap (x/d* = 0) is a 2-D wall value, so the 1-D edge is interpolated there.
    pts = {xd: p for xd, p in zip(d["xd"], d["pw"]) if p is not None and xd != 0.0}
    cp = 3.5 * AIR.R
    for st, row in d["surveys"]["adiabatic"].items():
        u = row[0]
        T = d["T0"] - u * u / (2 * cp)
        M = u / math.sqrt(1.4 * AIR.R * T)
        pts[st] = float(AIR.pe_pc_from_mach(M))
    xs = np.array(sorted(pts))
    lp = np.interp(x / d["dstar"], xs, np.log([pts[k] for k in xs]))
    M = np.sqrt(np.maximum(5.0 * (np.exp(-lp) ** (0.4 / 1.4) - 1.0), 1e-8))
    return Contour(x=x, r=r, theta=th, M=M), Rt


def _d5377_bl(inlet: str):
    """Two-layer BL from the measured nozzle-entrance state (TN D-5377 Table I(b))."""
    from rockmoc.analysis.boundary_layer import correct_wall_bl, two_layer_thicknesses

    d = _D5377
    c, Rt = _d5377_contour()
    xd = c.x / d["dstar"]
    Tw = np.interp(xd, d["xd"], d["tw_ad" if inlet == "adiabatic" else "tw_co"])
    th0, d0, D0 = d["inlet"][inlet]
    pt0 = two_layer_thicknesses(float(c.M[0]), Tw[0] / d["T0"], D0 / d0, 1.4)[1]
    bl = correct_wall_bl(c, AIR, throat_radius=Rt, pc=d["P0"], Tc=d["T0"], prandtl=0.71, wall_temperature=Tw,
                         viscosity=_sutherland_air, theta0=th0, phi0=th0 * pt0)
    return bl, xd


def _boundary_layer() -> list[Check]:
    d = _D5377
    out = []
    for inlet in ("adiabatic", "cooled"):
        bl, xd = _d5377_bl(inlet)
        for st, (_, t1, t2, s1, s2) in d["surveys"][inlet].items():
            i = int(np.argmin(np.abs(xd - st)))
            case = f"TN D-5377 {inlet} inlet, x/d*={st}"
            out.append(Check("conical (BL)", case, "θ [cm]", 0.5 * (t1 + t2), 100 * bl.theta_mom[i],
                             0.5 * abs(t2 - t1) / (0.5 * (t1 + t2)) + 0.20, "TN-D-5377 Tables II–VII",
                             note=f"band {t1}–{t2}; momentum+energy integral from inlet"))
            tiny = inlet == "cooled" and st < 1.0  # δ* ≈ 0 here: sign-sensitive, relative error meaningless
            out.append(Check("conical (BL)", case, "δ* [cm]", 0.5 * (s1 + s2), 100 * bl.delta_star[i], 0.35,
                             "TN-D-5377 Tables II–VII", informational=tiny,
                             note=("δ* ≈ 0 (thick cold thermal layer), relative error not meaningful; "
                                   if tiny else "") + f"band {min(s1, s2)}–{max(s1, s2)}"))
    return out


def _throat() -> list[Check]:
    from rockmoc.moc.transonic import HALL_MIN_RC_RT, hall_uv, sauer_offset, sauer_uv

    g = 1.4

    def p_of(V: float) -> float:
        return (1 - (g - 1) / (g + 1) * V * V) ** (g / (g - 1))

    def M_of(V: float) -> float:
        return math.sqrt(2 * V * V / ((g + 1) - (g - 1) * V * V))

    d = _D5377
    Rt, Rc = 0.5 * d["dstar"], d["Rc"]
    out = []
    for xd, pm in ((-0.117, .5733), (0.0, .4600)):  # Table I wall taps, Rc/Rt = 2.01
        x = xd * d["dstar"]
        rw = Rt + Rc - math.sqrt(Rc * Rc - x * x)
        out.append(Check("throat (Hall)", f"TN D-5377 wall, Rc/Rt=2.01, x/d*={xd}", "p_w/p0",
                         pm, p_of(math.hypot(*hall_uv(x, rw, g, Rt, Rc))), 0.02, "TN-D-5377 Table I"))
        us, vs = sauer_uv(x - sauer_offset(g, Rt, Rc), rw, g, Rt, Rc)
        out.append(Check("throat (Sauer)", f"TN D-5377 wall, Rc/Rt=2.01, x/d*={xd}", "p_w/p0",
                         pm, p_of(math.hypot(1 + us, vs)), 0.08, "TN-D-5377 Table I",
                         note="1st order; Hall preferred here", informational=True))
    # Cuffel: Rc/Rt = 0.625 < HALL_MIN_RC_RT, so the library uses Kliegel–Levine there
    from rockmoc.moc.transonic import discharge_coefficient, kl_uv

    for r, Mm, q, tol in ((0.0, 0.8, "axis", 0.05), (1.0, 1.4, "wall", 0.10)):
        out.append(Check("throat (KL)", f"Cuffel 1969 throat plane {q}, Rc/Rt=0.625", "M", Mm,
                         M_of(math.hypot(*kl_uv(0.0, r, g, 1.0, 0.625))), tol, "CUFFEL69 text",
                         note="Kliegel–Levine 3rd order"))
    out.append(Check("throat (KL)", "Cuffel 1969 flow coefficient, Rc/Rt=0.625", "Cd", 0.985,
                     discharge_coefficient(g, 0.625), 0.005, "CUFFEL69 text",
                     note="inviscid KL series; measured 0.985"))
    R = 0.625
    off = sauer_offset(g, 1.0, R)
    for r, Mm, q in ((0.0, 0.8, "axis"), (1.0, 1.4, "wall")):
        u, v = sauer_uv(-off, r, g, 1.0, R)
        out.append(Check("throat (Sauer)", f"Cuffel 1969 throat plane {q}, Rc/Rt=0.625", "M",
                         Mm, M_of(math.hypot(1 + u, v)), 0.30, "CUFFEL69 text", informational=True,
                         note=f"Sauer, for comparison (Hall diverges below Rc/Rt={HALL_MIN_RC_RT})"))
    return out


def _cone() -> list[Check]:
    out = []
    # Table I wall pressure vs 1-D isentropic at the local area ratio (d / d*)²
    for xd, dd, pm in ((3.665, 4.322, .0094), (5.497, 5.792, .0038), (7.394, 7.320, .0024), (0.818, 2.042, .1125)):
        eps = (dd / 1.492) ** 2
        M = float(AIR.mach_from_area_ratio(eps, branch="supersonic"))
        info = xd < 1.0
        out.append(Check("conical", f"TN D-5377 wall, x/d*={xd} (ε={eps:.1f})", "p_w/p0", pm,
                         float(AIR.pe_pc_from_mach(M)), 0.20, "TN-D-5377 Table I",
                         note="2-D near the arc/cone tangency" if info else "1-D", informational=info))
    return out


_TM107285_X = [0.000, 0.157, 0.175, 0.203, 0.231, 0.262, 0.295, 0.330, 0.366, 0.404, 0.444, 0.488, 0.533, 0.579,
               0.632, 0.686, 0.744, 0.805, 0.871, 0.940, 1.016, 1.095, 1.179, 1.267, 1.359, 1.450, 1.544, 1.623,
               1.727, 1.910, 2.217, 2.621, 2.891, 3.152, 3.495, 3.843, 4.280, 4.740, 5.314, 5.928, 6.302, 6.891,
               7.534, 8.037, 8.434, 8.989]
_TM107285_R = [0.638, 0.681, 0.691, 0.709, 0.726, 0.744, 0.765, 0.785, 0.808, 0.833, 0.856, 0.884, 0.912, 0.942,
               0.973, 1.006, 1.041, 1.080, 1.120, 1.161, 1.204, 1.252, 1.298, 1.349, 1.400, 1.450, 1.501, 1.542,
               1.595, 1.687, 1.836, 2.017, 2.129, 2.233, 2.362, 2.487, 2.631, 2.774, 2.936, 3.096, 3.185, 3.315,
               3.444, 3.538, 3.607, 3.693]
# SP-8120 chart, 80 % bell: ε → (θ_N, θ_e) [deg]; ε = 70 is extrapolated on the chart
_SP8120_80 = {10: (26.0, 10.0), 20: (28.5, 9.0), 40: (30.5, 8.0), 70: (33.0, 7.0)}


def _rao() -> list[Check]:
    from rockmoc.designers.rao_control import design_rao_control_surface
    from rockmoc.moc.rao import rao_length_reference

    out = []
    gas = PerfectGas(1.2, 350.0)
    for eps, (tn, te) in _SP8120_80.items():
        Me = float(gas.mach_from_area_ratio(eps, branch="supersonic"))
        _, m, _ = design_rao_control_surface(gas, throat_radius=0.05, exit_mach=Me, length_fraction=0.8, n_chars=80)
        note = "chart extrapolated" if eps > 50 else "chart read ±1°"
        out.append(Check("Rao TOC (MOC)", f"SP-8120 80 % bell, ε={eps}, γ=1.2", "θ_N [deg]", tn, m["theta_n_deg"],
                         1.5, "SP-8120 chart", rel=False, note=note))
        out.append(Check("Rao TOC (MOC)", f"SP-8120 80 % bell, ε={eps}, γ=1.2", "θ_e [deg]", te, m["theta_e_deg"],
                         1.5, "SP-8120 chart", rel=False, note=note))
    X, Rw = np.array(_TM107285_X), np.array(_TM107285_R)
    Rt = Rw[0]
    eps = (Rw[-1] / Rt) ** 2
    f = X[-1] / rao_length_reference(Rt, eps, 0.382 * Rt)
    Me = float(gas.mach_from_area_ratio(eps, branch="supersonic"))
    c, m, _ = design_rao_control_surface(gas, throat_radius=Rt, exit_mach=Me, length_fraction=f, n_chars=80)
    dev = float(np.max(np.abs(np.interp(X, c.x - c.x[0], c.r) - Rw)) / Rw[-1])
    th_e_ref = math.degrees(math.atan((Rw[-1] - Rw[-2]) / (X[-1] - X[-2])))
    out.append(Check("Rao TOC (MOC)", "TM-107285 Rao-code bell, ε=33.5, γ=1.2", "max |Δr| / r_e",
                     0.0, dev, 0.01, "TM-107285 Table I", rel=False, note="whole wall, 46 points"))
    out.append(Check("Rao TOC (MOC)", "TM-107285 Rao-code bell, ε=33.5, γ=1.2", "θ_e [deg]",
                     th_e_ref, m["theta_e_deg"], 1.0, "TM-107285 Table I", rel=False, note="last-segment slope"))
    return out


def _aerospike() -> list[Check]:
    from rockmoc.aerospike import design_spike

    out = []
    for g, eps in ((1.2, 10.0), (1.2, 18.0), (1.4, 8.0)):
        sp = design_spike(Pc=5e6, Tc=3300.0, gamma=g, molar_mass=22.0, Re=1.0, expansion_ratio=eps, n_points=4000)
        gas = PerfectGas(g)
        Me = float(gas.mach_from_area_ratio(eps, branch="supersonic"))
        nu_e = float(gas.prandtl_meyer(Me))
        Ms = np.linspace(1.0001, Me, 2000)
        mu = np.arcsin(1.0 / Ms)
        al = nu_e - gas.prandtl_meyer(Ms) + mu  # Mach-line angle from the lip
        # mass through the conical Mach-line surface: π(re²−r²) sinμ / sinα = (A/A*) A_t
        rr = np.sqrt(np.maximum(1.0 - gas.area_ratio(Ms) * np.sin(al) / (eps * np.sin(mu)), 0.0))
        xx = (1.0 - rr) / np.tan(al)
        m = (sp.X >= xx.min()) & (sp.X <= xx.max())
        dev = float(np.max(np.abs(sp.R[m] - np.interp(sp.X[m], xx, rr))))
        out.append(Check("aerospike", f"centred-fan plug, γ={g}, ε={eps}", "L / r_e", float(xx.max()),
                         float(sp.X.max()), 1e-3, "ANGELINO"))
        out.append(Check("aerospike", f"centred-fan plug, γ={g}, ε={eps}", "max |Δr| / r_e", 0.0, dev,
                         1e-4, "ANGELINO", rel=False))
    return out


# Zebbiche, Doctorat d'État thesis (Univ. Blida), Tables 5.42–5.47, perfect gas γ = 1.402,
# straight sonic line, lengths / throat radius: M_e → (L/y*, θ* corner, θ_max, y_e/y*)
_ZEBBICHE = {1.5: (2.0143, 2.4138, 3.7584, 1.0843), 2.0: (3.5731, 5.7466, 8.2830, 1.2984),
             3.0: (8.3554, 11.7091, 15.4652, 2.0542), 4.0: (16.8144, 16.1020, 20.2781, 3.2629),
             5.0: (30.7030, 19.2722, 23.5528, 4.9748), 6.0: (51.9449, 21.6141, 25.8751, 7.2440)}


def _mln_invariants() -> list[Check]:
    from rockmoc.designers import DesignSpec, design_mln

    out = []
    gas = PerfectGas(1.402, 287.05)
    for Me, (L, tc, tm, ye) in _ZEBBICHE.items():
        r = design_mln(DesignSpec(gas=gas, throat_radius=1.0, exit_mach=Me, n_characteristics=60,
                                  throat_mach=1.01), store_net=False)
        th = np.degrees(r.contour.theta)
        case = f"Zebbiche axisym MLN, γ=1.402, Me={Me}"
        src = "ZEBBICHE Tables 5.42–5.47"
        out += [Check("MLN (axisym)", case, "L / r*", L, r.metrics.length, 0.01, src),
                Check("MLN (axisym)", case, "corner θ* [deg]", tc, float(th[0]), 0.3, src, rel=False),
                Check("MLN (axisym)", case, "max wall θ [deg]", tm, float(th.max()), 0.3, src, rel=False,
                      note="max is downstream of the corner (inflection)"),
                Check("MLN (axisym)", case, "r_e / r*", ye, float(r.contour.r[-1]), 0.005, src)]
    for g, Me in ((1.4, 2.0), (1.4, 3.0), (1.2, 3.5)):
        gas = PerfectGas(g, 287.05)
        r = design_mln(DesignSpec(gas=gas, throat_radius=1.0, exit_mach=Me, n_characteristics=80), store_net=False)
        eps = float(gas.area_ratio(Me))
        out.append(Check("MLN (axisym)", f"γ={g}, Me={Me}", "A_e/A* (mass)", eps, r.metrics.area_ratio, 0.005,
                         "isentropic A/A*", note="mass-conservation consistency; no public axisym MLN table reachable"))
    return out


def _separation() -> list[Check]:
    """Kalt–Badal separated Cf against RocketCEA's independent implementation."""
    try:
        from rocketcea import separated_Cf as rc_sep
    except ImportError:
        return []
    from rockmoc.analysis.separation import separated_cf

    out = []
    for g, eps, pc, pa in ((1.25, 20.0, 200.0, 14.7), (1.2, 40.0, 1000.0, 14.7), (1.15, 60.0, 500.0, 14.7)):
        ref = rc_sep.ambientCf(gam=g, epsTot=eps, Pc=pc, Pamb=pa)[0]
        r = separated_cf(PerfectGas(g), area_ratio=eps, pc=pc, pa=pa, criterion="kalt_badal", recovery="kalt_badal")
        out.append(Check("separation", f"Kalt–Badal, γ={g}, ε={eps}, p0/pa={pc / pa:.1f}", "separated Cf", ref, r.Cf,
                         1e-4, "RocketCEA separated_Cf", note=f"ε_sep {r.eps_sep:.2f}; code-to-code"))
    return out


def _internal_shock() -> list[Check]:
    """Characteristic-crossing detector: TOP parabolas carry an internal shock, ideal walls do not
    (Frey et al., EUCASS 2017); conical nozzles with a circular-arc throat form one that reaches the
    centreline at z/r_th = 3.94 for Rc/Rt = 0.625 (Cuffel et al. 1969)."""
    from rockmoc.designers import DesignSpec, design_ideal, design_rao
    from rockmoc.moc.analysis import analyze_wall

    gas = PerfectGas(1.2, 350.0)
    spec = DesignSpec(gas=gas, throat_radius=1.0, exit_mach=3.5, n_characteristics=40)
    top = analyze_wall(design_rao(spec).contour, gas, rc_over_rt=1.5)
    ideal = analyze_wall(design_ideal(spec).contour, gas, rc_over_rt=1.5)
    Rc, a = 0.625, math.radians(15.0)
    xt = Rc * math.sin(a)
    x = np.linspace(0.0, 8.0, 1200)
    r = np.where(x <= xt, 1.0 + Rc - np.sqrt(np.maximum(Rc**2 - x**2, 0.0)),
                 1.0 + Rc * (1 - math.cos(a)) + (x - xt) * math.tan(a))
    cone = analyze_wall(Contour(x=x, r=r), AIR)
    src = "Frey et al. EUCASS 2017"
    out = [Check("internal shock", "TOP parabola, ε=13.8, γ=1.2", "shock detected", 1.0, float(top.has_internal_shock),
                 0.0, src, rel=False, note=f"origin x/Rt={top.shock_origin[0]:.2f}" if top.shock_origin else ""),
           Check("internal shock", "ideal contour, ε=13.8, γ=1.2", "shock detected", 0.0, float(ideal.has_internal_shock),
                 0.0, src, rel=False, note="shock-free by design")]
    if cone.shock_origin:
        out.append(Check("internal shock", "Cuffel 45–15° cone, Rc/Rt=0.625", "first crossing x / r_th", 3.94,
                         cone.shock_origin[0], 0.20, "CUFFEL69 text", informational=True,
                         note="reference is where the shock meets the centreline; crossing found near the axis"))
    return out


def _kinetics() -> list[Check]:
    """Finite-rate nozzle: equilibrium end vs NASA CEA, and the frozen / fast-rate limits."""
    try:
        from rocketcea.cea_obj import CEA_Obj

        from rockmoc.thermo.kinetics import kinetic_nozzle
    except ImportError:
        return []
    from rockmoc.designers import DesignSpec, design_rao_control

    c = design_rao_control(DesignSpec(gas=PerfectGas(1.2, 520.0), throat_radius=1.0, exit_mach=4.2,
                                      n_characteristics=40)).contour
    x, r = (c.x - c.x[0]) * 0.05, c.r * 0.05
    eps = float((r[-1] / r.min()) ** 2)
    k = kinetic_nozzle(reactants="H2:3.0, O2:1.0", pc=3e6, x=x, r=r)
    k0 = kinetic_nozzle(reactants="H2:3.0, O2:1.0", pc=3e6, x=x, r=r, rate_multiplier=0.0)
    kf = kinetic_nozzle(reactants="H2:3.0, O2:1.0", pc=3e6, x=x, r=r, rate_multiplier=1e3)
    cea = CEA_Obj(oxName="GOX", fuelName="GH2").get_Isp(Pc=3e6 / 6894.757, MR=32.0 / 6.0, eps=eps, frozen=0)
    case = f"GH2/GO2 MR 5.33, pc 3 MPa, ε={eps:.1f}, Rt 5 cm"
    return [
        Check("kinetics", case, "shifting-eq Isp_vac [s]", cea, k.Isp_shifting, 0.005, "NASA CEA (RocketCEA)"),
        Check("kinetics", case, "zero-rate Isp = frozen [s]", k0.Isp_frozen, k0.Isp_kinetic, 1e-4, "frozen limit",
              note="checks the derived quasi-1D equations"),
        Check("kinetics", case, "fast-rate Isp → shifting [s]", kf.Isp_shifting, kf.Isp_kinetic, 2e-3,
              "equilibrium limit", note=f"η_kin at real rates {k.eta_kinetic:.4f}"),
    ]


def _two_phase() -> list[Check]:
    """Gas–particle lag model: tiny particles must recover the closed-form mixture isentrope."""
    from rockmoc.designers import DesignSpec, design_rao_control
    from rockmoc.performance.two_phase import two_phase_nozzle

    c = design_rao_control(DesignSpec(gas=PerfectGas(1.2, 350.0), throat_radius=1.0, exit_mach=3.0,
                                      n_characteristics=40)).contour
    x, r = (c.x - c.x[0]) * 0.05, c.r * 0.05
    kw = dict(x=x, r=r, gamma_gas=1.2, R_gas=350.0, pc=7e6, Tc=3400.0, particle_mass_fraction=0.3)
    nolag = two_phase_nozzle(diameter=1e-8, **kw)
    typ = two_phase_nozzle(diameter=4e-6, **kw)
    return [
        Check("two-phase", "30 % Al2O3, d → 0, Rt 5 cm", "Isp_lag / Isp_no-lag", 1.0, nolag.eta_two_phase, 1e-4,
              "equilibrium mixture isentrope", note="checks the derived gas–particle equations"),
        Check("two-phase", "30 % Al2O3, d = 4 µm (≈1 µm per inch of throat), Rt 5 cm", "η_2φ", 0.975,
              typ.eta_two_phase, 0.02, "typical 2–4 % loss (Sutton; US pat. 8051640 rule)", informational=True,
              note="plausibility only; Hermsen/SPP correlations not verified"),
    ]


def _euler_throat() -> list[Check]:
    """2-D Euler throat solver (Richardson-extrapolated) on Cuffel's Rc/Rt = 0.625 nozzle."""
    from rockmoc.moc.transonic_euler import cuffel_wall, solve_throat_extrapolated

    r = solve_throat_extrapolated(cuffel_wall(0.625, r_inlet=2.0), -1.8, 1.0)
    case, src = "Cuffel 1969 throat plane, Rc/Rt=0.625", "CUFFEL69 text"
    return [Check("throat (Euler)", case, "axis M", 0.8, r["M_axis"], 0.05, src),
            Check("throat (Euler)", case, "wall M", 1.4, r["M_wall"], 0.10, src),
            Check("throat (Euler)", case, "Cd", 0.985, r["Cd"], 0.005, src, note="inviscid; measured 0.985")]


def run_experimental_validation() -> list[Check]:
    """All literature / experiment comparisons (≈1 min)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MOCFallbackWarning)
        return [*_boundary_layer(), *_throat(), *_cone(), *_rao(), *_aerospike(), *_mln_invariants(),
                *_separation(), *_internal_shock(), *_kinetics(), *_two_phase(), *_euler_throat()]


def experimental_report_markdown(rows: list[Check]) -> str:
    gated = [c for c in rows if not c.informational]
    lines = [
        "# Validation against published data",
        "",
        f"Generated by `python -m rockmoc.validate.experimental`. "
        f"Gated checks: {sum(c.passed for c in gated)}/{len(gated)} pass; "
        f"{len(rows) - len(gated)} informational rows (known model limits).",
        "",
        "| Nozzle | Case | Quantity | Reference | rockmoc | Error | Tol | Pass | Source | Note |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for c in rows:
        tol = f"±{100 * c.tol:.0f} %" if c.rel else f"±{c.tol:g}"
        flag = "info" if c.informational else ("✅" if c.passed else "❌")
        lines.append(f"| {c.nozzle} | {c.case} | {c.quantity} | {c.reference:.4g} | {c.predicted:.4g} | "
                     f"{c.error} | {tol} | {flag} | {c.source} | {c.note} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    from pathlib import Path

    rows = run_experimental_validation()
    md = experimental_report_markdown(rows)
    out = Path(__file__).resolve().parents[3] / "docs" / "validation_experimental.md"
    out.write_text(md, encoding="utf-8")
    print(md)
