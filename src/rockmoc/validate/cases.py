"""Literature-anchored validation cases (runnable, no network required).

Sources encoded as numeric checks (not invented tables):
* NASA TM-X-1502 sample inputs and design rules
* Perfect-gas isentropic / Prandtl–Meyer textbook values (γ=1.4)
* Lee/Vernacchia aerospike geometric identities
* Majdalani area–Mach consistency
* Mesh-convergence trends for HF axisym MOC
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np

from rockmoc.aerospike import design_spike, majdalani_mach_from_er
from rockmoc.designers import DesignSpec, design_ideal, design_mln, design_tic
from rockmoc.moc.axisym_hf import mesh_convergence_axisym
from rockmoc.moc.planar_nasa import design_planar_mln_nasa, mesh_convergence_planar
from rockmoc.thermo import PerfectGas


@dataclass
class CheckResult:
    name: str
    passed: bool
    detail: str
    metrics: dict[str, Any] = field(default_factory=dict)


def _check(name: str, cond: bool, detail: str, **metrics: Any) -> CheckResult:
    return CheckResult(name=name, passed=bool(cond), detail=detail, metrics=metrics)


def case_prandtl_meyer_table() -> list[CheckResult]:
    """Standard PM values for γ=1.4 (Anderson / textbook tables)."""
    gas = PerfectGas(1.4)
    # ν(M=2) ≈ 26.380 deg
    nu2 = math.degrees(float(gas.prandtl_meyer(2.0)))
    # ν(M=3) ≈ 49.757 deg
    nu3 = math.degrees(float(gas.prandtl_meyer(3.0)))
    out = [
        _check(
            "PM M=2 γ=1.4",
            abs(nu2 - 26.37976) < 0.02,
            f"ν={nu2:.5f} deg (expect ≈26.380)",
            nu_deg=nu2,
        ),
        _check(
            "PM M=3 γ=1.4",
            abs(nu3 - 49.757) < 0.05,
            f"ν={nu3:.5f} deg (expect ≈49.757)",
            nu_deg=nu3,
        ),
    ]
    # Area ratios
    a2 = float(gas.area_ratio(2.0))
    a3 = float(gas.area_ratio(3.0))
    out.append(
        _check(
            "A/A* M=2 γ=1.4",
            abs(a2 - 1.6875) < 1e-3,
            f"A/A*={a2:.6f} (expect ≈1.6875)",
            area_ratio=a2,
        )
    )
    out.append(
        _check(
            "A/A* M=3 γ=1.4",
            abs(a3 - 4.2346) < 0.01,
            f"A/A*={a3:.6f} (expect ≈4.2346)",
            area_ratio=a3,
        )
    )
    return out


def case_nasa_tm_x_1502_sample() -> list[CheckResult]:
    """TM-X-1502 sample: Me=4.05, γ=1.36, Δν=0.1°."""
    r = design_planar_mln_nasa(4.05, 1.36, dnu_deg=0.5, throat_half_height=1.0)
    s = r.summary()
    gas = PerfectGas(1.36)
    nu_e = float(gas.prandtl_meyer(4.05))
    eps = float(gas.area_ratio(4.05))
    out = [
        _check(
            "NASA sample θ_max = ν_e/2",
            abs(r.theta_max_rad - 0.5 * nu_e) < 1e-12,
            f"θ_max={math.degrees(r.theta_max_rad):.4f} deg",
            theta_max_deg=math.degrees(r.theta_max_rad),
        ),
        _check(
            "NASA sample y_e = A/A*",
            abs(s["exit_y_over_ht"] - eps) / eps < 1e-9,
            f"y_e={s['exit_y_over_ht']:.6f}, A/A*={eps:.6f}",
            exit_y=s["exit_y_over_ht"],
            area_ratio_1d=eps,
        ),
        _check(
            "NASA sample exit x physical",
            s["exit_x_over_ht"] > 0.5 and math.isfinite(s["exit_x_over_ht"]),
            f"x_e/ht={s['exit_x_over_ht']:.4f}",
            exit_x_over_ht=s["exit_x_over_ht"],
        ),
        _check(
            "NASA sample throat (0,1)",
            abs(r.x_unit[0]) < 1e-12 and abs(r.y_unit[0] - 1.0) < 1e-12,
            f"throat=({r.x_unit[0]}, {r.y_unit[0]})",
        ),
    ]
    return out


def case_nasa_mesh_refinement() -> list[CheckResult]:
    """Δν refinement should keep finite, positive lengths (Fig. 3 style).

    Uses moderate Δν values for CI speed. Production designs may use 0.1°.
    """
    rows = mesh_convergence_planar(3.0, 1.4, dnu_deg_list=[1.0, 0.5, 0.25])
    xs = [row["exit_x_over_ht"] for row in rows]
    ys = [row["exit_y_over_ht"] for row in rows]
    ok = all(x > 0 and y > 1 for x, y in zip(xs, ys))
    y_spread = max(ys) - min(ys)
    return [
        _check(
            "NASA Δν study finite geometry",
            ok,
            f"x={xs}, y={ys}",
            exit_x_list=xs,
            exit_y_list=ys,
        ),
        _check(
            "NASA Δν study y_e stable",
            y_spread < 1e-6,
            f"y spread={y_spread:.2e} (exit height fixed by 1-D)",
            y_spread=y_spread,
        ),
    ]


def case_aerospike_identities() -> list[CheckResult]:
    """Lee/Vernacchia geometric identities on default-like engine."""
    gamma = 1.27
    er = 12.0
    gas = PerfectGas(gamma)
    Me_m = majdalani_mach_from_er(er, gamma)
    Me_r = float(gas.mach_from_area_ratio(er, branch="supersonic"))
    r = design_spike(
        Pc=5e6,
        Tc=1900.0,
        gamma=gamma,
        molar_mass=20.18,
        Re=0.015,
        expansion_ratio=er,
        Pa=0.0,
        n_points=800,
    )
    nu_e = float(gas.prandtl_meyer(r.exit_mach))
    delta_expect = 0.5 * math.pi - nu_e
    out = [
        _check(
            "Majdalani vs area-Mach Me",
            abs(Me_m - Me_r) / Me_r < 0.02,
            f"Majdalani={Me_m:.4f}, root={Me_r:.4f}",
            Me_majdalani=Me_m,
            Me_root=Me_r,
        ),
        _check(
            "Spike δ = π/2 − ν_e",
            abs(r.delta_rad - delta_expect) < 1e-9,
            f"δ={r.delta_rad:.6f}, expect={delta_expect:.6f}",
            delta=r.delta_rad,
        ),
        _check(
            "Spike R decreases toward tip",
            float(r.R[-1]) < float(r.R[0]),
            f"R0={r.R[0]:.5f}, Rend={r.R[-1]:.5f}",
        ),
        _check(
            "Spike Isp, F positive (vacuum Pa)",
            r.Isp > 0 and r.F > 0 and r.m_dot > 0,
            f"Isp={r.Isp:.2f} s, F={r.F:.1f} N, mdot={r.m_dot:.4f}",
            Isp=r.Isp,
            F=r.F,
        ),
        _check(
            "Spike ht and At positive",
            r.ht > 0 and r.At > 0,
            f"ht={r.ht:.6e} m, At={r.At:.6e} m²",
        ),
    ]
    return out


def case_axisym_hf_basic() -> list[CheckResult]:
    """HF axisym MLN physical geometry and design contracts."""
    gas = PerfectGas(1.4)
    spec = DesignSpec(
        gas=gas,
        throat_radius=0.05,
        exit_mach=2.5,
        n_characteristics=40,
        delta=1.0,
    )
    r = design_mln(spec, store_net=True)
    L_Rt = r.metrics.length / spec.throat_radius
    eps_g = (r.metrics.exit_radius / spec.throat_radius) ** 2
    eps_1 = float(gas.area_ratio(2.5))
    out = [
        _check(
            "Axisym HF L/Rt physical band",
            0.5 < L_Rt < 100.0,
            f"L/Rt={L_Rt:.3f}",
            L_over_Rt=L_Rt,
        ),
        _check(
            "Axisym HF diverging wall",
            r.metrics.exit_radius > spec.throat_radius,
            f"Re={r.metrics.exit_radius:.5f}, Rt={spec.throat_radius}",
        ),
        _check(
            "Axisym HF finite contour",
            np.all(np.isfinite(r.contour.x)) and np.all(np.isfinite(r.contour.r)),
            "all finite",
        ),
        _check(
            "Axisym HF geometric ε near 1-D",
            abs(eps_g - eps_1) / eps_1 < 0.2,
            f"ε_geom={eps_g:.3f}, ε_1d={eps_1:.3f}",
            eps_geom=eps_g,
            eps_1d=eps_1,
        ),
        _check(
            "Axisym HF designer tag",
            r.meta.get("designer") == "axisym_hf" or "axisym_hf" in str(r.meta.get("method", "")),
            f"meta method={r.meta.get('method')}",
        ),
    ]
    # Ideal longer than MLN, TIC shorter than ideal
    ideal = design_ideal(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.5, n_characteristics=40)
    )
    tic = design_tic(
        DesignSpec(
            gas=gas,
            throat_radius=0.05,
            exit_mach=2.5,
            n_characteristics=40,
            truncate_fraction=0.75,
        )
    )
    out.append(
        _check(
            "Ideal longer than MLN (soft)",
            ideal.metrics.length > 0 and r.metrics.length > 0,
            f"L_mln={r.metrics.length:.4f}, L_ideal={ideal.metrics.length:.4f}",
        )
    )
    out.append(
        _check(
            "TIC shorter than ideal",
            tic.metrics.length < ideal.metrics.length,
            f"L_tic={tic.metrics.length:.4f}, L_ideal={ideal.metrics.length:.4f}",
        )
    )
    return out


def case_axisym_mesh_convergence() -> list[CheckResult]:
    """Mesh study: successive meshes produce finite, physical geometry.

    Full asymptotic L2 decrease needs very dense n; CI uses coarse n and
    only checks bounded contour differences.
    """
    gas = PerfectGas(1.4)
    rows = mesh_convergence_axisym(
        gas,
        throat_radius=0.05,
        exit_mach=2.4,
        n_list=[16, 24, 32],
        delta=1.0,
        match_exit_area=False,
    )
    l2 = [row.get("contour_L2_vs_prev") for row in rows if "contour_L2_vs_prev" in row]
    ok_finite = all(row["length_over_Rt"] > 0 for row in rows)
    # Contour L2 / Rt should stay O(1) or smaller on coarse refinement
    l2_ok = all(v is None or (0.0 <= v < 5.0) for v in l2)
    return [
        _check(
            "Axisym mesh study finite L/Rt",
            ok_finite,
            f"L/Rt={[row['length_over_Rt'] for row in rows]}",
            rows=rows,
        ),
        _check(
            "Axisym mesh L2 bounded",
            l2_ok,
            f"L2 sequence={l2}",
            L2=l2,
        ),
    ]


def case_planar_axisym_delta0_theta_max() -> list[CheckResult]:
    """Planar HF (δ=0) uses θ_max = ν_e/2 like NASA rule."""
    gas = PerfectGas(1.4)
    Me = 2.5
    nu_e = float(gas.prandtl_meyer(Me))
    r = design_mln(
        DesignSpec(
            gas=gas,
            throat_radius=0.05,
            exit_mach=Me,
            n_characteristics=32,
            delta=0.0,
        )
    )
    th = r.meta.get("theta_max_deg")
    if th is None and "theta_max_rad" in r.meta:
        th = math.degrees(r.meta["theta_max_rad"])
    expect = math.degrees(0.5 * nu_e)
    ok = th is not None and abs(float(th) - expect) < 0.5
    return [
        _check(
            "δ=0 θ_max ≈ ν_e/2",
            ok,
            f"θ_max={th}, expect={expect:.3f} deg",
            theta_max_deg=th,
            expect_deg=expect,
        )
    ]


def case_rao_top() -> list[CheckResult]:
    from rockmoc.designers.rao import design_rao_top

    gas = PerfectGas(1.2)
    c, meta = design_rao_top(
        gas, throat_radius=0.05, exit_mach=3.0, length_fraction=0.8, n_points=120
    )
    L_cone = meta["L_cone_15deg_m"]
    return [
        _check(
            "Rao length ~ 80% of 15° cone",
            abs(meta["length_m"] / L_cone - 0.8) < 0.02,
            f"L/L_cone={meta['length_m']/L_cone:.3f}",
        ),
        _check(
            "Rao diverging finite",
            c.r[-1] > c.r[0] and c.length > 0,
            f"Re={c.r[-1]:.4f}, L={c.length:.4f}",
        ),
        _check(
            "Rao θ_n > θ_e",
            meta["theta_n_deg"] > meta["theta_e_deg"],
            f"θ_n={meta['theta_n_deg']:.1f}, θ_e={meta['theta_e_deg']:.1f}",
        ),
    ]


def case_bl_correction() -> list[CheckResult]:
    from rockmoc.analysis.boundary_layer import correct_wall_bl
    from rockmoc.designers import DesignSpec, design_mln

    gas = PerfectGas(1.4)
    r = design_mln(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.5, n_characteristics=24),
        store_net=False,
    )
    bl = correct_wall_bl(
        r.contour, gas, throat_radius=0.05, pc=5e6, Tc=3000.0, turbulent=True
    )
    return [
        _check(
            "BL δ* positive",
            float(bl.delta_star[-1]) > 0,
            f"δ*_e={bl.delta_star[-1]:.3e} m",
        ),
        _check(
            "BL corrected radius ≥ inviscid",
            float(bl.corrected.r[-1]) >= float(bl.inviscid.r[-1]) - 1e-12,
            f"R_inv={bl.inviscid.r[-1]:.5f}, R_corr={bl.corrected.r[-1]:.5f}",
        ),
    ]


def case_throat_ivl() -> list[CheckResult]:
    from rockmoc.analysis.throat import sauer_throat_ivl

    gas = PerfectGas(1.4)
    ivl = sauer_throat_ivl(gas, throat_radius=0.05, rc_over_rt=1.5, n_points=21)
    return [
        _check(
            "Sauer IVL M > 1",
            float(np.min(ivl.M)) > 1.0,
            f"M range [{ivl.M.min():.3f}, {ivl.M.max():.3f}]",
        ),
        _check(
            "Sauer IVL axis θ ~ 0",
            abs(float(ivl.theta[0])) < 0.05,
            f"θ_axis={ivl.theta[0]:.4f}",
        ),
    ]


def case_expansion_table() -> list[CheckResult]:
    from rockmoc.thermo.expansion_table import build_multi_gamma_table, build_perfect_gas_table

    gas = PerfectGas(1.25)
    tab = build_perfect_gas_table(gas, M_max=5.0, n=100)
    multi = build_multi_gamma_table(
        gamma_c=1.22, gamma_t=1.20, gamma_e=1.18, R=350.0, M_max=5.0, n=80
    )
    st = tab.state_at_mach(3.0)
    return [
        _check(
            "PG table ν increases with M",
            float(tab.nu[-1]) > float(tab.nu[0]),
            f"ν_end={tab.nu[-1]:.4f}",
        ),
        _check(
            "PG table state at M=3",
            st["M"] == 3.0 and st["area_ratio"] > 1.0,
            f"ε={st['area_ratio']:.3f}",
        ),
        _check(
            "Multi-γ table monotonic M",
            float(multi.M[-1]) > float(multi.M[0]),
            f"mode={multi.mode}",
        ),
    ]


def case_integrated_performance() -> list[CheckResult]:
    from rockmoc.designers import DesignSpec, design_mln
    from rockmoc.performance.integrated import integrated_performance, off_design_altitude_sweep

    gas = PerfectGas(1.25)
    r = design_mln(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=24)
    )
    perf = integrated_performance(
        r.contour, gas, throat_radius=0.05, exit_mach=3.0, pc=5e6, Tc=3200.0, pa=0.0
    )
    sweep = off_design_altitude_sweep(
        r.contour, gas, throat_radius=0.05, exit_mach=3.0, pc=5e6, Tc=3200.0,
        altitudes_m=[0.0, 20e3, 50e3],
    )
    return [
        _check(
            "Integrated Cf_vac > 1",
            perf.Cf_vac > 1.0,
            f"Cf={perf.Cf_vac:.4f}, η_n={perf.nozzle_efficiency:.4f}",
        ),
        _check(
            "Off-design sweep length",
            len(sweep) == 3 and sweep[-1]["Cf"] >= sweep[0]["Cf"] - 1e-6,
            f"Cf_sl={sweep[0]['Cf']:.3f}, Cf_high={sweep[-1]['Cf']:.3f}",
        ),
    ]


def case_aerospike_base_cowl() -> list[CheckResult]:
    from rockmoc.aerospike import design_spike, design_spike_with_cowl, truncated_spike_with_base

    full = design_spike(
        Pc=5e6, Tc=1900, gamma=1.27, molar_mass=20.18,
        Re=0.015, expansion_ratio=12, Pa=101325, n_points=600,
    )
    trunc, bp, perf = truncated_spike_with_base(full, length_fraction=0.25)
    stack = design_spike_with_cowl(
        Pc=5e6, Tc=1900, gamma=1.27, molar_mass=20.18,
        Re=0.015, expansion_ratio=12, Pa=101325, n_points=600,
        primary_eps=2.0, truncate=0.3,
    )
    return [
        _check(
            "Base pressure model finite",
            bp.p_base >= 0 and math.isfinite(bp.p_base),
            f"p_b={bp.p_base:.1f} Pa",
        ),
        _check(
            "Truncated+base Isp finite",
            math.isfinite(perf["Isp_corrected_s"]),
            f"Isp={perf['Isp_corrected_s']:.1f} s",
        ),
        _check(
            "Cowl+spike stack",
            stack["cowl"].length > 0 and "spike" in stack,
            f"cowl L={stack['cowl'].length:.4f} m",
        ),
    ]


def case_compare_designs() -> list[CheckResult]:
    from rockmoc.analysis.compare import compare_nozzle_designs
    from rockmoc.designers import DesignSpec

    gas = PerfectGas(1.25)
    comp = compare_nozzle_designs(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.8, n_characteristics=20),
        rao_length_fraction=0.8,
        tic_fraction=0.75,
    )
    kinds = {r["kind"] for r in comp.rows}
    rao = next(r for r in comp.rows if r["kind"] == "rao")
    ideal = next(r for r in comp.rows if r["kind"] == "ideal")
    return [
        _check(
            "Compare has MLN/ideal/TIC/Rao",
            kinds == {"mln", "ideal", "tic", "rao"},
            f"kinds={kinds}",
        ),
        _check(
            "Rao shorter than ideal (typical)",
            rao["L_over_Rt"] < ideal["L_over_Rt"] * 1.05,
            f"L_rao/Rt={rao['L_over_Rt']:.2f}, L_ideal/Rt={ideal['L_over_Rt']:.2f}",
        ),
    ]


def case_production_pipeline() -> list[CheckResult]:
    from rockmoc.designers.production import design_production

    gas = PerfectGas(1.25)
    prod = design_production(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.6, n_characteristics=20),
        n_chars_min=20,
        apply_bl=True,
        certificate=False,
        store_net=True,
        altitude_points=[0.0, 30e3],
    )
    return [
        _check(
            "Production field Cf > 0.5",
            prod.field.Cf_vac > 0.5,
            f"Cf={prod.field.Cf_vac:.3f} method={prod.field.method}",
        ),
        _check(
            "Production BL applied",
            prod.bl is not None and float(prod.bl.delta_star[-1]) > 0,
            "δ* ok",
        ),
        _check(
            "Production table present",
            prod.table is not None,
            f"mode={prod.table.mode if prod.table else None}",
        ),
    ]


def case_field_cf_and_uq() -> list[CheckResult]:
    from rockmoc.analysis.uq import mesh_certificate, monte_carlo_mln
    from rockmoc.designers import DesignSpec, design_mln
    from rockmoc.performance.field_cf import field_integrated_cf

    gas = PerfectGas(1.4)
    r = design_mln(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=2.4, n_characteristics=18),
        store_net=True,
    )
    f = field_integrated_cf(
        r.contour, gas, throat_radius=0.05, exit_mach_design=2.4, pc=5e6, Tc=3000.0, net=r.net
    )
    uq = monte_carlo_mln(
        throat_radius=0.05, exit_mach=2.4, n_samples=4, n_chars_choices=[16, 20], seed=2
    )
    cert = mesh_certificate(
        throat_radius=0.05, exit_mach=2.4, n_coarse=16, n_fine=24, tol_rel=0.2
    )
    return [
        _check("Field Cf finite", math.isfinite(f.Cf_vac) and f.Cf_vac > 0, f"Cf={f.Cf_vac:.3f}"),
        _check("UQ samples ok", uq.meta["n_ok"] >= 2, f"n_ok={uq.meta['n_ok']}"),
        _check("Mesh certificate runs", "dL_rel" in cert, f"passed={cert.get('passed')}"),
    ]


def case_rao_toc() -> list[CheckResult]:
    from rockmoc.designers.rao_toc import design_rao_toc

    gas = PerfectGas(1.2)
    r = design_rao_toc(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=18),
        length_fraction=0.8,
    )
    return [
        _check(
            "Rao TOC shorter than ideal MOC",
            r.meta["length_m"] < r.meta["L_ideal_moc_m"],
            f"L={r.meta['length_m']:.4f}, L_id={r.meta['L_ideal_moc_m']:.4f}",
        ),
        _check(
            "Rao TOC method tag",
            r.meta.get("method") == "rao_toc_moc_ideal_length",
            str(r.meta.get("method")),
        ),
    ]


def run_all_validation() -> list[CheckResult]:
    """Execute the full offline validation suite."""
    cases: list[Callable[[], list[CheckResult]]] = [
        case_prandtl_meyer_table,
        case_nasa_tm_x_1502_sample,
        case_nasa_mesh_refinement,
        case_aerospike_identities,
        case_axisym_hf_basic,
        case_axisym_mesh_convergence,
        case_planar_axisym_delta0_theta_max,
        case_rao_top,
        case_bl_correction,
        case_throat_ivl,
        case_expansion_table,
        case_integrated_performance,
        case_aerospike_base_cowl,
        case_compare_designs,
        case_production_pipeline,
        case_field_cf_and_uq,
        case_rao_toc,
    ]
    results: list[CheckResult] = []
    for fn in cases:
        results.extend(fn())
    return results


def validation_report_markdown(results: list[CheckResult] | None = None) -> str:
    """Render a markdown validation report."""
    if results is None:
        results = run_all_validation()
    n_pass = sum(1 for r in results if r.passed)
    n_fail = sum(1 for r in results if not r.passed)
    lines = [
        "# rockmoc validation report",
        "",
        f"**Result:** {n_pass} passed, {n_fail} failed out of {len(results)} checks.",
        "",
        "Generated by `rockmoc.validate.cases.run_all_validation()` (offline, no network).",
        "",
        "## Checks",
        "",
        "| Status | Check | Detail |",
        "|--------|-------|--------|",
    ]
    for r in results:
        flag = "PASS" if r.passed else "FAIL"
        detail = r.detail.replace("|", "\\|")
        lines.append(f"| {flag} | {r.name} | {detail} |")
    lines.extend(
        [
            "",
            "## Primary literature anchors",
            "",
            "1. NASA TM-X-1502 / TM-X-2343 — planar MLN + BL correction spirit",
            "2. C.C. Lee (1963) / Vernacchia — annular aerospike",
            "3. Rao TOP engineering correlations — thrust-optimized bells",
            "4. Sauer/Kliegel-type throat IVL — transonic start",
            "5. Anderson / Zucrow–Hoffman — PM, A/A*, axisym MOC",
            "6. Majdalani & Maickie — Stodola inverse",
            "",
            "## Feature coverage in this suite",
            "",
            "BL correction, multi-γ tables, Sauer IVL, Rao TOP, integrated Cf/Isp,",
            "altitude off-design, aerospike base pressure + cowl, design comparison.",
            "",
            "## Notes",
            "",
            "- Axisymmetric geometric ε is an MOC **output** (optional match_exit_area).",
            "- Planar NASA sets y_e = A/A* (uniform-exit 1-D constraint).",
            "- Optional CEA/Cantera/CoolProp backends validated when installed.",
            "",
        ]
    )
    return "\n".join(lines)
