"""MOC-aware / contour-integrated performance and off-design analysis."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from rockmoc.geometry.contour import Contour
from rockmoc.thermo.perfect_gas import PerfectGas

G0 = 9.80665


@dataclass(slots=True)
class IntegratedPerformance:
    """Thrust coefficient and efficiency estimates from contour + 1-D edge."""

    Cf_vac: float
    Cf_amb: float
    Isp_vac: float | None
    Isp_amb: float | None
    pe_pc: float
    Me: float
    area_ratio: float
    divergence_efficiency: float
    kinetic_efficiency: float
    nozzle_efficiency: float
    meta: dict

    def summary(self) -> dict:
        return {
            "Cf_vac": self.Cf_vac,
            "Cf_amb": self.Cf_amb,
            "Isp_vac_s": self.Isp_vac,
            "Isp_amb_s": self.Isp_amb,
            "pe_pc": self.pe_pc,
            "Me": self.Me,
            "area_ratio": self.area_ratio,
            "divergence_efficiency": self.divergence_efficiency,
            "kinetic_efficiency": self.kinetic_efficiency,
            "nozzle_efficiency": self.nozzle_efficiency,
            **self.meta,
        }


def divergence_efficiency(contour: Contour) -> float:
    """Exit-plane divergence efficiency from exit wall angle.

    Standard conical nozzle correction (also used for contoured bells when
    only θ_e is known)::

        η_div = (1 + cos θ_e) / 2

    For a well-designed MLN / ideal bell, θ_e ≈ 0 ⇒ η_div ≈ 1.
    """
    if contour.theta is None or len(contour.theta) < 1:
        return 1.0
    th_e = abs(float(contour.theta[-1]))
    return float(0.5 * (1.0 + math.cos(th_e)))


def exit_kinetic_efficiency(contour: Contour, *, n_quad: int = 32) -> float:
    """Approximate axial kinetic efficiency at the exit plane.

    Assumes a linear θ(r) from 0 on axis to θ_e at the wall and uniform |V|,
    then η_ke = (∫ u² r dr) / (∫ V² r dr) with u = V cos θ.

    For θ_e = 0 this is exactly 1 (no artificial wall-history penalty).
    """
    if contour.theta is None or len(contour.theta) < 1:
        return 1.0
    th_e = abs(float(contour.theta[-1]))
    if th_e < 1e-12:
        return 1.0
    Re = float(contour.r[-1])
    if Re <= 0.0:
        return 1.0
    # ∫_0^{Re} cos²(θ_e r/Re) r dr  /  ∫ r dr
    r = np.linspace(0.0, Re, n_quad)
    th = th_e * (r / Re)
    trapz = getattr(np, "trapezoid", None) or np.trapz
    num = float(trapz(np.cos(th) ** 2 * r, r))
    den = float(trapz(r, r))
    if den <= 0.0:
        return 1.0
    return float(np.clip(num / den, 0.90, 1.0))


_MOC_CACHE: dict = {}


def moc_thrust_efficiency(
    contour: Contour, gas, *, rc_over_rt: float | None = None, n_chars: int = 60
) -> dict:
    """Inviscid geometric thrust efficiency of a *given* wall from a direct MOC march.

    Momentum balance on start line → wall → exit plane (axisymmetric, vacuum)::

        F = ∫_start (ρu_x² + p) dA + ∫_wall p dA_x,     ṁ = ∫_start ρ u_x dA

    using ρu² = γpM² and ρu = γpM/a, so only M, θ along the start line and the wall
    pressure are needed. Wall stretch not reached by the march (the last C+ leave
    through the exit) gets the 1-D pressure scaled to the last marched wall point.
    η = (F/ṁ) / (F/ṁ)_1-D at the same ε. Captures over-expansion and lip
    recompression of short bells that a θ_e-only loss misses (CFD: Rao 60 %).
    Also returns the wall p/p0(ε) used by the separation model.
    """
    key = (contour.x.tobytes(), contour.r.tobytes(), id(gas), n_chars, rc_over_rt)
    if key in _MOC_CACHE:
        return _MOC_CACHE[key]
    from rockmoc.moc.analysis import analyze_wall

    wa = analyze_wall(contour, gas, n_chars=n_chars, rc_over_rt=rc_over_rt)
    g0 = float(gas.gamma)

    def st(M: float) -> tuple[float, float, float]:  # p/p0, γ_local, a/a0
        g = float(gas.gamma_at_mach(M))
        return float(gas.pe_pc_from_mach(M)), g, math.sqrt(g / g0 / float(gas.temperature_ratio(M)))

    s0 = min(a + b for a, b in wa.points)
    line = sorted((p for (a, b), p in wa.points.items() if a + b == s0), key=lambda p: p.r)
    F = md = 0.0
    for p1, p2 in zip(line[:-1], line[1:], strict=True):
        f = []
        for q in (p1, p2):
            pp, g, a = st(q.M)
            c = math.cos(q.theta)
            f.append((g * pp * q.M**2 * c * c + pp, g * pp * q.M * c / a))
        dA = math.pi * (p2.r**2 - p1.r**2)
        F += 0.5 * (f[0][0] + f[1][0]) * dA
        md += 0.5 * (f[0][1] + f[1][1]) * dA
    i_t = int(np.argmin(contour.r))
    Rt, x_t = float(contour.r[i_t]), float(contour.x[i_t])
    rr = [line[-1].r] + [q.r for q in wa.wall]
    pw = [st(line[-1].M)[0]] + [st(q.M)[0] for q in wa.wall]
    Mw = [line[-1].M] + [q.M for q in wa.wall]
    x_last = wa.wall[-1].x + x_t if wa.wall else float(contour.x[i_t])

    def p1d(r: float) -> float:
        return float(
            gas.pe_pc_from_mach(
                gas.mach_from_area_ratio(max((r / Rt) ** 2, 1.0 + 1e-9), branch="supersonic")
            )
        )

    tail = contour.x > x_last + 1e-12 * max(abs(x_last), 1.0)
    if tail.any():
        k = pw[-1] / p1d(rr[-1])
        for r_ in contour.r[tail]:
            rr.append(float(r_))
            pw.append(k * p1d(float(r_)))
            Mw.append(
                float(
                    gas.mach_from_area_ratio(max((r_ / Rt) ** 2, 1.0 + 1e-9), branch="supersonic")
                )
            )
    rr_a, pw_a = np.array(rr), np.array(pw)
    F += float(np.sum(0.5 * (pw_a[1:] + pw_a[:-1]) * math.pi * np.diff(rr_a**2)))
    eps = (float(contour.r[-1]) / Rt) ** 2
    Me = float(gas.mach_from_area_ratio(eps, branch="supersonic"))
    _, g, a = st(Me)
    ideal = a * (g * Me**2 + 1.0) / (g * Me)
    out = {
        "eta": (F / md) / ideal,
        "coverage": float(wa.meta["wall_coverage"]),
        "has_internal_shock": wa.has_internal_shock,
        "wall_eps": (rr_a / Rt) ** 2,
        "wall_p_pc": pw_a,
        "wall_M": np.array(Mw),
        "n_chars": n_chars,
    }
    if len(_MOC_CACHE) > 64:
        _MOC_CACHE.clear()
    _MOC_CACHE[key] = out
    return out


def _is_straight_exit(contour: Contour, tol_deg: float = 0.2) -> bool:
    """Cone-like wall: θ constant over the last half of the length (source flow, λ exact)."""
    if contour.theta is None or len(contour.theta) < 6:
        return False
    x = np.asarray(contour.x)
    th = np.degrees(np.asarray(contour.theta)[x >= x[0] + 0.5 * (x[-1] - x[0])])
    return bool(th.size >= 3 and np.ptp(th) < tol_deg)


def integrated_performance(
    contour: Contour,
    gas: PerfectGas,
    *,
    throat_radius: float,
    exit_mach: float,
    pc: float = 5e6,
    Tc: float = 3000.0,
    pa: float = 0.0,
    delta_cf_friction: float = 0.0,
    eta_chemical: float = 1.0,
    eta_two_phase: float = 1.0,
    exit_displacement_area: float = 0.0,
    separation: str | None = None,
    separation_recovery: str = "kalt_badal",
    geometric: str = "auto",
    rc_over_rt: float | None = None,
) -> IntegratedPerformance:
    """Combine 1-D isentropic Cf with a JANNAF-style nozzle loss stack.

    ``Cf = Cf_1d · η_geom · η_chem · η_2φ − ΔCf_friction``

    η_geom is ONE geometric (divergence) factor, chosen by ``geometric``:

    * ``"auto"`` (default): θ_e < 0.5° → 1 (designer exit is uniform); straight
      (conical) exit → η_div = (1+cos θ_e)/2, exact for source flow; any other
      wall (TOP/TOC bells, truncated contours) → :func:`moc_thrust_efficiency`
      (direct MOC march on the wall), falling back to η_div if the march fails.
    * ``"divergence"``: η_div only.  ``"moc"``: force the MOC march.

    v0.9 multiplied η_div · η_ke, which counts the same exit-angle loss twice
    (CFD: cones and bells under-predicted by 1.4–3.8 % Isp). η_ke is still
    reported in ``loss_stack`` but not applied. Wall-angle history is
    *not* used as a loss proxy. ``delta_cf_friction`` comes from
    :attr:`BLCorrectionResult.delta_cf_friction`; ``eta_chemical`` (finite-rate
    / frozen-vs-equilibrium, see
    :func:`rockmoc.thermo.rocketcea_backend.kinetic_efficiency_bounds`)
    and ``eta_two_phase`` are user inputs (1.0 = no loss).

    ``separation`` (a criterion name from :mod:`rockmoc.analysis.separation`,
    e.g. ``"schmucker"``) replaces the attached-flow ambient Cf by the
    free-shock-separated value when the nozzle is over-expanded enough to separate.
    """
    for name, eta in (("eta_chemical", eta_chemical), ("eta_two_phase", eta_two_phase)):
        if not 0.0 < eta <= 1.0:
            raise ValueError(f"{name} must be in (0, 1]; got {eta}")
    if delta_cf_friction < 0.0:
        raise ValueError("delta_cf_friction must be ≥ 0")
    Rt = throat_radius
    Re = float(contour.r[-1])
    eps = (Re / Rt) ** 2
    Me = float(exit_mach)
    pe_pc = float(gas.pe_pc_from_mach(Me))
    pa_pc = pa / pc if pc > 0 else 0.0

    Cf_1d_vac = gas.thrust_coefficient(pe_pc, 0.0, exit_mach=Me, area_ratio=eps)

    if geometric not in ("auto", "divergence", "moc"):
        raise ValueError(f"geometric must be 'auto', 'divergence' or 'moc'; got {geometric!r}")
    eta_div = divergence_efficiency(contour)
    eta_ke = exit_kinetic_efficiency(contour)
    th_e = (
        abs(float(contour.theta[-1])) if contour.theta is not None and len(contour.theta) else 0.0
    )
    moc = None
    if geometric == "divergence":
        eta_geom, geom_model = eta_div, "divergence"
    elif geometric == "auto" and th_e < math.radians(0.5):
        eta_geom, geom_model = 1.0, "uniform_exit"
    elif geometric == "auto" and _is_straight_exit(contour):
        eta_geom, geom_model = eta_div, "conical_source_flow"
    else:
        try:
            moc = moc_thrust_efficiency(contour, gas, rc_over_rt=rc_over_rt)
            if moc["coverage"] < 0.6 or not 0.85 < moc["eta"] <= 1.003:
                raise ValueError(
                    f"MOC march unusable (coverage {moc['coverage']:.2f}, eta {moc['eta']:.4f})"
                )
            eta_geom, geom_model = min(moc["eta"], 1.0), "moc_wall_integral"
        except Exception as exc:  # noqa: BLE001 - any march failure → closed-form fallback
            if geometric == "moc":
                raise
            moc, eta_geom, geom_model = None, eta_div, f"divergence (MOC fallback: {exc})"
    eta_mult = eta_geom * eta_chemical * eta_two_phase

    Cf_vac = Cf_1d_vac * eta_mult - delta_cf_friction
    Cf_amb = Cf_vac - pa_pc * (eps + exit_displacement_area / (math.pi * throat_radius**2))
    eta_n = Cf_vac / Cf_1d_vac if Cf_1d_vac > 0 else eta_mult
    sep_meta = None
    if separation and pa > 0.0:
        from rockmoc.analysis.separation import separated_cf

        wall = (
            {"wall_eps": moc["wall_eps"], "wall_p_pc": moc["wall_p_pc"], "wall_M": moc["wall_M"]}
            if moc
            else {}
        )
        sep = separated_cf(
            gas,
            area_ratio=eps,
            pc=pc,
            pa=pa,
            criterion=separation,  # type: ignore[arg-type]
            recovery=separation_recovery,
            **wall,
        )  # type: ignore[arg-type]
        sep_meta = {
            "criterion": separation,
            "separated": sep.separated,
            "eps_sep": sep.eps_sep,
            "psep_pa": sep.psep_pa,
            "full_flow_tp1207": sep.full_flow_tp1207,
            "recovery": separation_recovery,
            "wall_pressure": "moc" if moc else "1-D",
        }
        if sep.separated:
            # losses act on the attached (vacuum) part only; the ambient/recovery term is a pressure force
            Cf_amb = (
                sep.Cf_vac_attached * eta_mult + (sep.Cf - sep.Cf_vac_attached) - delta_cf_friction
            )

    cstar = gas.characteristic_velocity(Tc)
    Isp_vac = Cf_vac * cstar / G0
    Isp_amb = Cf_amb * cstar / G0

    return IntegratedPerformance(
        Cf_vac=float(Cf_vac),
        Cf_amb=float(Cf_amb),
        Isp_vac=float(Isp_vac),
        Isp_amb=float(Isp_amb),
        pe_pc=pe_pc,
        Me=Me,
        area_ratio=float(eps),
        divergence_efficiency=eta_div,
        kinetic_efficiency=eta_ke,
        nozzle_efficiency=eta_n,
        meta={
            "Cf_1d_vac": float(Cf_1d_vac),
            "cstar_m_s": float(cstar),
            "pc_Pa": pc,
            "Tc_K": Tc,
            "pa_Pa": pa,
            "separation": sep_meta,
            "geometric_model": geom_model,
            "moc_internal_shock": bool(moc["has_internal_shock"]) if moc else None,
            "loss_stack": {
                "eta_geometric": eta_geom,
                "eta_divergence": eta_div,
                "eta_exit_kinetic_energy": eta_ke,
                "eta_chemical": eta_chemical,
                "eta_two_phase": eta_two_phase,
                "eta_friction": 1.0 - delta_cf_friction / (Cf_1d_vac * eta_mult)
                if Cf_1d_vac > 0
                else 1.0,
                "delta_cf_friction": delta_cf_friction,
            },
            "theta_exit_deg": float(
                math.degrees(abs(float(contour.theta[-1]))) if contour.theta is not None else 0.0
            ),
        },
    )


def off_design_altitude_sweep(
    contour: Contour,
    gas: PerfectGas,
    *,
    throat_radius: float,
    exit_mach: float,
    pc: float = 5e6,
    Tc: float = 3000.0,
    altitudes_m: list[float] | None = None,
    **losses: float | str | None,
) -> list[dict]:
    """Cf / Isp vs altitude using simple atmosphere Pa(h).

    ``losses`` are forwarded to :func:`integrated_performance`
    (``delta_cf_friction``, ``eta_chemical``, ``eta_two_phase``, ``separation``).
    """
    from rockmoc.analysis.atmosphere import pressure_at_altitude

    if altitudes_m is None:
        altitudes_m = [0.0, 5e3, 10e3, 15e3, 20e3, 30e3, 50e3, 80e3]
    rows = []
    for h in altitudes_m:
        pa = pressure_at_altitude(h)
        perf = integrated_performance(
            contour, gas, throat_radius=throat_radius, exit_mach=exit_mach, pc=pc, Tc=Tc, pa=pa,
            **losses,
        )
        rows.append(
            {
                "altitude_m": h,
                "pa_Pa": pa,
                "Cf": perf.Cf_amb,
                "Isp_s": perf.Isp_amb,
                "pe_pc": perf.pe_pc,
                "pa_pc": pa / pc,
            }
        )
    return rows
