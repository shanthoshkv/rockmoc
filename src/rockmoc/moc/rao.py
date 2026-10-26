"""Rao thrust-optimised contour (TOC) by the method of characteristics.

Construction (Rao 1958; Zucrow & Hoffman Vol. 2 §16.6), vacuum optimum:

1. Smooth-throat expansion: circular arc of radius Rc up to the inflection
   angle θ_N, solved by the same characteristic lattice as the ideal contour.
2. ND is the last C− from the inflection point N to the axis.
3. Control surface DE is a C+ from a point D on ND to the exit lip E, along
   which Rao's optimality conditions hold::

       V cos(θ − μ) / cos μ        = λ   (const)
       r^δ ρ V² sin²θ tan μ        = K   (const)

   and the mass flow across DE equals the mass flow across DN.
4. D is chosen so the lip radius r_E equals the design exit radius; θ_N so
   the lip lies at the requested length.
5. The wall NE is the streamline through N in the region between ND and DE
   (Goursat problem + mass-flow wall, as for the ideal contour).

Velocities are normalised by a* (V = M*) and densities by ρ*, so ρV = 1/(A/A*).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.optimize import brentq

from rockmoc.exceptions import DesignError, MeshError, ThermoError, warn_fallback
from rockmoc.geometry.contour import Contour
from rockmoc.moc._common import (
    _fp,
    _poly_mass,
    _records_fallbacks,
    _wall_to_contour,
    blend_contours,
)
from rockmoc.moc.ideal import build_ideal_nets
from rockmoc.moc.transition import _transition_wall
from rockmoc.moc.transonic import mstar_from_mach
from rockmoc.moc.types import CharNet, FlowPoint, PointKind
from rockmoc.thermo.gas_protocol import MOCGas


def _lerp_point(a: FlowPoint, b: FlowPoint, f: float, gas: MOCGas) -> FlowPoint:
    return _fp(a.x + f * (b.x - a.x), a.r + f * (b.r - a.r), a.theta + f * (b.theta - a.theta),
               a.nu + f * (b.nu - a.nu), gas, PointKind.INTERIOR)


def _control_surface(
    D: FlowPoint, m_target: float, gas: MOCGas, *, Rt: float, delta: float, n_steps: int = 400
) -> list[FlowPoint] | None:
    """March Rao's control surface (a C+) from D until it carries ``m_target``.

    Returns the points of DE (D first, lip E last) or None if the Rao
    conditions cannot be satisfied (e.g. D on the axis, where θ = 0).
    """
    g_of = gas.gamma_at_mach  # local γ for M*
    planar = abs(delta) < 0.5
    if D.theta <= 1e-6 or D.r <= 0.0:
        return None

    def _res(th: float, M: float, r: float) -> tuple[float, float]:
        # λ- and K-conditions as smooth functions of (θ, M). Eliminating θ via
        # θ = μ ± acos(λ cos μ / M*) is singular at θ = μ, which the optimum
        # DE crosses for L/L15 ≈ 0.8 bells; the 2×2 Jacobian is regular there.
        mu = math.asin(1.0 / M)
        ms = mstar_from_mach(M, float(g_of(M)))
        e1 = ms * math.cos(th - mu) / math.cos(mu) - lam
        e2 = (r if not planar else 1.0) * ms * ms / float(gas.area_ratio(M)) \
            * math.sin(th) ** 2 * math.tan(mu) / K - 1.0
        return e1, e2

    def _solve(th: float, M: float, r: float) -> tuple[float, float] | None:
        for _ in range(40):
            f1, f2 = _res(th, M, r)
            if abs(f1) < 1e-13 and abs(f2) < 1e-13:
                return th, M
            ht, hm = 1e-7, 1e-7 * M
            a1, a2 = _res(th + ht, M, r)
            b1, b2 = _res(th, M + hm, r)
            j11, j21 = (a1 - f1) / ht, (a2 - f2) / ht
            j12, j22 = (b1 - f1) / hm, (b2 - f2) / hm
            det = j11 * j22 - j12 * j21
            if det == 0.0 or not math.isfinite(det):
                return None
            dth = (f1 * j22 - f2 * j12) / det
            dM = (j11 * f2 - j21 * f1) / det
            lim = 1.0
            while M - lim * dM <= 1.0 + 1e-9 or not 0.0 < th - lim * dth < 0.5 * math.pi:
                lim *= 0.5
                if lim < 1e-6:
                    return None
            th, M = th - lim * dth, M - lim * dM
            if abs(dth) < 1e-14 and abs(dM) < 1e-14 * M:
                return th, M
        return None

    ms_d = mstar_from_mach(D.M, float(g_of(D.M)))
    lam = ms_d * math.cos(D.theta - D.mu) / math.cos(D.mu)
    K = (D.r if not planar else 1.0) * ms_d * ms_d / float(gas.area_ratio(D.M)) \
        * math.sin(D.theta) ** 2 * math.tan(D.mu)

    pts = [D]
    m_acc, x, r, M = 0.0, D.x, D.r, D.M
    th, mu = D.theta, D.mu
    ar = float(gas.area_ratio(M))
    dr = Rt / n_steps * 4.0
    dth_prev = dM_prev = 0.0
    for _ in range(50 * n_steps):
        r_new = r + dr
        # Newton on (θ, M) from a linear extrapolation of the last step
        st = _solve(min(max(th + dth_prev, 1e-6), 1.5), max(M + dM_prev, 1.0 + 1e-6), r_new)
        if st is None:
            return None
        th_n, M_new = st
        dth_prev, dM_prev = th_n - th, M_new - M
        mu_n = math.asin(1.0 / M_new)
        ar_n = float(gas.area_ratio(M_new))
        phi = 0.5 * (th + mu + th_n + mu_n)
        if math.sin(phi) <= 0.0 or math.tan(phi) <= 0.0:
            return None
        x_new = x + dr / math.tan(phi)
        w = 1.0 if planar else math.pi * (r + r_new)
        dm = w * 0.5 * (math.sin(mu) / ar + math.sin(mu_n) / ar_n) * dr / math.sin(phi)
        if m_acc + dm >= m_target:
            f = (m_target - m_acc) / dm
            x_e, r_e = x + f * (x_new - x), r + f * dr
            st_e = _solve(th + f * (th_n - th), M + f * (M_new - M), r_e)
            if st_e is None:
                return None
            pts.append(_fp(x_e, r_e, st_e[0], float(gas.prandtl_meyer(st_e[1])), gas, PointKind.INTERIOR))
            return pts
        m_acc += dm
        x, r, M, th, mu, ar = x_new, r_new, M_new, th_n, mu_n, ar_n
        pts.append(_fp(x, r, th, float(gas.prandtl_meyer(M)), gas, PointKind.INTERIOR))
    return None


def _thin(pts: list[FlowPoint], n: int) -> list[FlowPoint]:
    if len(pts) <= n + 1:
        return pts
    idx = np.unique(np.linspace(0, len(pts) - 1, n + 1).round().astype(int))
    return [pts[i] for i in idx]


def _lattice(gas: MOCGas, *, Rt: float, Rc: float, theta_n: float, n_chars: int,
             throat_mach: float, delta: float) -> tuple[dict, str]:
    for start in ("hall", "sauer", "uniform"):
        lat: dict = {}
        try:
            build_ideal_nets(gas, Rt=Rt, Rc=Rc, theta_attach=theta_n, n_chars=n_chars,
                             throat_mach=throat_mach, delta=delta, transition=False,
                             initial_line=start, lattice=lat)  # type: ignore[arg-type]
            return lat, start
        except (MeshError, ThermoError):
            continue
    raise DesignError(f"throat expansion to θ_N = {math.degrees(theta_n):.2f}° failed")


def _solve_d(ND: list[FlowPoint], mND: list[float], r_exit: float, gas: MOCGas, *,
             Rt: float, delta: float) -> tuple[float, list[FlowPoint]] | None:
    """Position t along ND (0 = axis, len−1 = N) with lip radius r_E = r_exit."""
    mdot = mND[-1]

    def _de(t: float) -> list[FlowPoint] | None:
        i = min(int(t), len(ND) - 2)
        f = t - i
        D = _lerp_point(ND[i], ND[i + 1], f, gas)
        m_d = mND[i] + f * (mND[i + 1] - mND[i])
        return _control_surface(D, mdot - m_d, gas, Rt=Rt, delta=delta)

    def _h(t: float) -> float:
        de = _de(t)
        return math.nan if de is None else de[-1].r - r_exit

    ts = np.linspace(0.5, len(ND) - 1.05, 24)
    hs = [_h(float(t)) for t in ts]
    for k in range(len(ts) - 1):
        a, b = hs[k], hs[k + 1]
        if not (math.isnan(a) or math.isnan(b)) and a * b <= 0.0:
            t = float(brentq(_h, float(ts[k]), float(ts[k + 1]), xtol=1e-9))
            de = _de(t)
            return (t, de) if de is not None else None
    return None


def rao_length_reference(Rt: float, eps: float, Rc: float) -> float:
    """15° conical nozzle length incl. throat arc (Rao's length-fraction basis)."""
    a15 = math.radians(15.0)
    return (Rt * (math.sqrt(eps) - 1.0) + Rc * (1.0 / math.cos(a15) - 1.0)) / math.tan(a15)


@_records_fallbacks
def design_rao_moc(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    length_fraction: float = 0.80,
    n_chars: int = 40,
    rc_over_rt: float = 0.382,
    throat_mach: float = 1.02,
    delta: float = 1.0,
    clamp_length: bool = False,
) -> tuple[Contour, CharNet, dict[str, Any]]:
    """Rao TOC for the 1-D exit area ratio of ``exit_mach`` and length L = f · L₁₅.

    If the length is not reachable on the lattice, raise — or with
    ``clamp_length`` use the nearest feasible length and record a warning.
    """
    if exit_mach <= 1.0 or throat_radius <= 0.0:
        raise DesignError("need exit_mach > 1 and throat_radius > 0")
    Rt = float(throat_radius)
    Rc = rc_over_rt * Rt
    planar = abs(delta) < 0.5
    eps = float(gas.area_ratio(float(exit_mach)))
    r_exit = Rt * (eps if planar else math.sqrt(eps))
    L_target = length_fraction * rao_length_reference(Rt, eps, Rc)
    nu_e = float(gas.prandtl_meyer(float(exit_mach)))

    # One lattice serves every candidate inflection node N: the flow upstream
    # of a node's C− does not depend on where the arc ends.
    lat, start = _lattice(gas, Rt=Rt, Rc=Rc, theta_n=min(0.5 * nu_e, math.radians(45.0)),
                          n_chars=n_chars, throat_mach=throat_mach, delta=delta)
    P, m, b_max = lat["P"], lat["m"], lat["b_end"]

    def _node(b: int) -> Any:
        if (b + m, b) not in P:
            return None
        ND = [P[(a, b)] for a in range(b + m, b - m - 1, -1)]
        mND = _poly_mass(ND, gas, planar=planar)
        sol = _solve_d(ND, mND, r_exit, gas, Rt=Rt, delta=delta)
        return None if sol is None else (b, ND, mND, sol)

    nodes = {b: s_ for b in range(m + 1, b_max + 1) if (s_ := _node(b)) is not None}
    L15 = rao_length_reference(Rt, eps, Rc)
    pick = None
    for b in sorted(nodes):
        if b + 1 not in nodes:  # only interpolate between adjacent arc nodes
            continue
        lo, hi = nodes[b], nodes[b + 1]
        g_lo, g_hi = lo[3][1][-1].x - L_target, hi[3][1][-1].x - L_target
        if g_lo * g_hi <= 0.0:
            pick = (lo, hi, g_lo / (g_lo - g_hi) if g_lo != g_hi else 1.0)
            break
    if pick is None and clamp_length:
        if nodes:
            best = min(nodes, key=lambda b: abs(nodes[b][3][1][-1].x - L_target))
            lo = nodes[best]
            hi = nodes.get(best + 1, lo)
            pick = (lo, hi, 0.0)
            warn_fallback(
                f"L/L15 = {length_fraction:.2f} not reachable; used "
                f"{lo[3][1][-1].x / L15:.3f} (nearest feasible)"
            )
    if pick is None:
        f_ok = sorted(s_[3][1][-1].x / L15 for s_ in nodes.values())
        span = f"{f_ok[0]:.2f}–{f_ok[-1]:.2f}" if f_ok else "none"
        raise DesignError(
            f"Rao TOC: L/L15 = {length_fraction:.2f} not reachable at ε = {eps:.2f} "
            f"(feasible on this lattice: {span}; try n_chars ↑ or another length)"
        )

    def _build(node: Any, c_minus: list | None = None) -> tuple[list[FlowPoint], list[list[FlowPoint]], list[FlowPoint]]:
        b, ND, mND, (t_d, DE) = node
        a_star = b - m
        arc = [P[(aw, aw + m)] for aw in range(0, a_star + 1)]
        # C+ rows start on ND at points spaced uniformly in mass between D and N
        # (ND nodes crowd near N, which would leave the wall under-resolved).
        i = min(int(t_d), len(ND) - 2)
        m_d = mND[i] + (t_d - i) * (mND[i + 1] - mND[i])
        bases = [ND[: i + 1] + [DE[0]]]
        for mk in np.linspace(m_d, mND[-1], max(2 * len(ND), 2 * n_chars))[1:-1]:
            k = max(int(np.searchsorted(mND, mk)) - 1, 0)
            f = (mk - mND[k]) / max(mND[k + 1] - mND[k], 1e-300)
            bases.append(ND[: k + 1] + [_lerp_point(ND[k], ND[k + 1], f, gas)])
        t_wall, t_rows = _transition_wall(gas, bases, Rt=Rt, delta=delta, mdot=mND[-1],
                                          exit_line=_thin(DE, max(len(bases), 16)),
                                          c_minus_out=c_minus)
        return arc + [w for w in t_wall if w.x > arc[-1].x], t_rows, DE

    lo, hi, w = pick
    wall_lo, _, de_lo = _build(lo)
    net = CharNet()
    wall, t_rows, DE = _build(hi, net.c_minus)
    b_end = hi[0]
    theta_n = (1.0 - w) * P[(lo[0] - m, lo[0])].theta + w * P[(b_end - m, b_end)].theta

    net.wall = wall
    net.axis = [P[(bb + m, bb)] for bb in range(0, b_end + 1) if (bb + m, bb) in P]
    rows: dict[int, list[tuple[int, FlowPoint]]] = {}
    for (a, b), p in P.items():
        if b <= b_end and a <= b_end + m:
            rows.setdefault(a, []).append((b, p))
    net.interior = [[p for _, p in sorted(v, key=lambda t: t[0])] for _, v in sorted(rows.items())] + t_rows
    cols: dict[int, list[tuple[int, FlowPoint]]] = {}
    for (a, b), p in P.items():
        if b <= b_end and a <= b_end + m:
            cols.setdefault(b, []).append((a, p))
    net.c_minus += [[p for _, p in sorted(v, key=lambda t: t[0])] for _, v in sorted(cols.items()) if len(v) >= 2]
    # Lattice points past the exit plane are outside the nozzle: trim them
    x_lip = wall[-1].x * (1.0 + 1e-9)
    inside = lambda ln: [p for p in ln if p.x <= x_lip]  # noqa: E731
    net.interior = [r_ for r_ in map(inside, net.interior) if len(r_) >= 2]
    net.c_minus = [c_ for c_ in map(inside, net.c_minus) if len(c_) >= 2]
    net.axis = inside(net.axis)
    net.points = [p for r_ in net.interior for p in r_]

    contour = blend_contours(_wall_to_contour(wall_lo, Rt, axial_exit=False),
                             _wall_to_contour(wall, Rt, axial_exit=False), w,
                             n=max(200, 6 * n_chars), axial_exit=False)
    E_th = (1.0 - w) * de_lo[-1].theta + w * DE[-1].theta
    E_M = (1.0 - w) * de_lo[-1].M + w * DE[-1].M
    meta = {
        "method": "rao_moc_toc",
        "theta_n_deg": math.degrees(theta_n),
        "theta_e_deg": math.degrees(E_th),
        "exit_lip_mach": E_M,
        "length_m": float(contour.length),
        "L15_m": rao_length_reference(Rt, eps, Rc),
        "length_fraction": length_fraction,
        "area_ratio_target": eps,
        "rc_over_rt": rc_over_rt,
        "initial_line": start,
        "reference": "Rao 1958 control-surface TOC via MOC",
    }
    return contour, net, meta
