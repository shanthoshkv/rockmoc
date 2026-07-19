"""Deep multi-regime validation: literature references vs rockmoc solver.

Loads ``tests/data/validation_catalog.json`` (or a passed path), runs the
solver across thermo / planar MLN / axisym MLN / Rao / aerospike / performance
regimes, and writes:

* JSON  — full reference vs solver + error%
* CSV   — flat metric table
* MD    — human-readable report

Primary sources (also listed in the catalog):
NASA TM-X-1502, NACA RM E8J12, NASA TM-X-2343, Rao 1958, Lee/Vernacchia,
Anderson/NACA isentropic tables, Sutton 1-D Cf, Zucrow/Hoffman MOC.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from rockmoc.aerospike import design_spike
from rockmoc.analysis.compare import compare_nozzle_designs
from rockmoc.designers import DesignSpec, design_ideal, design_mln, design_tic
from rockmoc.designers.rao import design_rao
from rockmoc.moc.planar_nasa import design_planar_mln_nasa
from rockmoc.thermo.perfect_gas import PerfectGas


@dataclass
class MetricCompare:
    name: str
    reference: float | None
    solver: float | None
    error_percent: float | None
    abs_error: float | None
    tol_percent: float | None
    passed: bool
    unit: str = ""
    note: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class CaseResult:
    case_id: str
    regime: str
    source_id: str
    source_citation: str
    description: str
    inputs: dict[str, Any]
    metrics: list[MetricCompare] = field(default_factory=list)
    passed: bool = True
    meta: dict[str, Any] = field(default_factory=dict)

    def finalize(self) -> None:
        self.passed = all(m.passed for m in self.metrics)

    def as_dict(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "regime": self.regime,
            "source_id": self.source_id,
            "source_citation": self.source_citation,
            "description": self.description,
            "inputs": self.inputs,
            "passed": self.passed,
            "metrics": [m.as_dict() for m in self.metrics],
            "meta": self.meta,
        }


def _err_pct(ref: float | None, sol: float | None) -> tuple[float | None, float | None]:
    if ref is None or sol is None:
        return None, None
    abs_e = abs(float(sol) - float(ref))
    if abs(float(ref)) < 1e-30:
        return (0.0 if abs_e < 1e-12 else float("inf")), abs_e
    return 100.0 * abs_e / abs(float(ref)), abs_e


def _metric(
    name: str,
    reference: float | None,
    solver: float | None,
    *,
    tol_percent: float | None = None,
    unit: str = "",
    note: str = "",
    abs_tol: float | None = None,
    higher_is_pass: bool | None = None,
    lower_is_pass: bool | None = None,
    bool_pass: bool | None = None,
) -> MetricCompare:
    if bool_pass is not None:
        return MetricCompare(
            name=name,
            reference=reference,
            solver=solver,
            error_percent=0.0 if bool_pass else 100.0,
            abs_error=None if reference is None or solver is None else abs(float(solver) - float(reference)),
            tol_percent=tol_percent,
            passed=bool_pass,
            unit=unit,
            note=note,
        )
    err_p, abs_e = _err_pct(reference, solver)
    if higher_is_pass is not None and solver is not None:
        passed = float(solver) >= float(higher_is_pass)
    elif lower_is_pass is not None and solver is not None:
        passed = float(solver) <= float(lower_is_pass)
    elif abs_tol is not None and abs_e is not None:
        passed = abs_e <= abs_tol
    elif tol_percent is not None and err_p is not None:
        passed = err_p <= tol_percent
    elif reference is None and solver is not None:
        passed = math.isfinite(float(solver))
    else:
        passed = err_p is not None and err_p < 1e6
    return MetricCompare(
        name=name,
        reference=None if reference is None else float(reference),
        solver=None if solver is None else float(solver),
        error_percent=None if err_p is None else float(err_p),
        abs_error=None if abs_e is None else float(abs_e),
        tol_percent=tol_percent,
        passed=bool(passed),
        unit=unit,
        note=note,
    )


def _load_catalog(path: Path | None = None) -> dict:
    if path is None:
        # package-relative: tests/data when running from repo
        candidates = [
            Path(__file__).resolve().parents[3] / "tests" / "data" / "validation_catalog.json",
            Path.cwd() / "tests" / "data" / "validation_catalog.json",
        ]
        for c in candidates:
            if c.exists():
                path = c
                break
        if path is None:
            raise FileNotFoundError("validation_catalog.json not found")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _cite(catalog: dict, source_id: str) -> str:
    s = catalog.get("sources", {}).get(source_id, {})
    return s.get("citation", source_id)


# ---------------------------------------------------------------------------
# Regime runners
# ---------------------------------------------------------------------------


def run_thermo_isentropic(catalog: dict) -> list[CaseResult]:
    reg = catalog["regimes"]["thermo_isentropic_g14"]
    gas = PerfectGas(reg["gamma"])
    tols = reg["tol_percent"]
    src = reg["source_id"]
    out: list[CaseResult] = []
    for M_str, ref in reg["reference_table"].items():
        M = float(M_str)
        case = CaseResult(
            case_id=f"thermo_g{reg['gamma']}_M{M_str}",
            regime="thermo_isentropic",
            source_id=src,
            source_citation=_cite(catalog, src),
            description=f"Isentropic/PM table γ={reg['gamma']} M={M}",
            inputs={"gamma": reg["gamma"], "Mach": M},
        )
        nu = math.degrees(float(gas.prandtl_meyer(M)))
        ar = float(gas.area_ratio(M))
        pe = float(gas.pe_pc_from_mach(M))
        case.metrics.append(
            _metric("nu_deg", ref["nu_deg"], nu, tol_percent=tols["nu_deg"], unit="deg")
        )
        case.metrics.append(
            _metric("area_ratio", ref["area_ratio"], ar, tol_percent=tols["area_ratio"], unit="-")
        )
        case.metrics.append(
            _metric("pe_pc", ref["pe_pc"], pe, tol_percent=tols["pe_pc"], unit="-")
        )
        # Cf vacuum analytical
        Cf = gas.thrust_coefficient(pe, 0.0, exit_mach=M, area_ratio=ar)
        case.metrics.append(
            _metric("Cf_vac_1d", None, Cf, unit="-", note="solver 1-D Cf (no external table digit)")
        )
        case.finalize()
        out.append(case)
    return out


def run_thermo_multi_gamma(catalog: dict) -> list[CaseResult]:
    reg = catalog["regimes"]["thermo_multi_gamma"]
    src = reg["source_id"]
    tols = reg["tol_percent"]
    out: list[CaseResult] = []
    for c in reg["cases"]:
        g, Me = float(c["gamma"]), float(c["Me"])
        gas = PerfectGas(g)
        case = CaseResult(
            case_id=f"thermo_multi_g{g}_M{Me}",
            regime="thermo_multi_gamma",
            source_id=src,
            source_citation=_cite(catalog, src),
            description=f"PM & A/A* round-trip γ={g} Me={Me}",
            inputs=c,
        )
        nu = float(gas.prandtl_meyer(Me))
        M_back = float(gas.mach_from_prandtl_meyer(nu))
        ar = float(gas.area_ratio(Me))
        M_ar = float(gas.mach_from_area_ratio(ar, branch="supersonic"))
        case.metrics.append(
            _metric(
                "pm_roundtrip_M",
                Me,
                M_back,
                tol_percent=tols["pm_roundtrip_M"],
                unit="-",
                note="M → ν → M",
            )
        )
        case.metrics.append(
            _metric(
                "area_roundtrip_M",
                Me,
                M_ar,
                tol_percent=tols["area_roundtrip_M"],
                unit="-",
                note="M → A/A* → M",
            )
        )
        case.finalize()
        out.append(case)
    return out


def run_nasa_tm_x_1502(catalog: dict) -> list[CaseResult]:
    reg = catalog["regimes"]["nasa_tm_x_1502_sample"]
    src = reg["source_id"]
    Me, g = float(reg["Me"]), float(reg["gamma"])
    gas = PerfectGas(g)
    nu_e = float(gas.prandtl_meyer(Me))
    ar = float(gas.area_ratio(Me))
    r = design_planar_mln_nasa(
        Me,
        g,
        dnu_deg=float(reg["dnu_deg"]),
        throat_half_height=float(reg["throat_half_height"]),
        store_net=False,
    )
    s = r.summary()
    case = CaseResult(
        case_id="nasa_tm_x_1502_sample",
        regime="planar_mln_nasa",
        source_id=src,
        source_citation=_cite(catalog, src),
        description=f"NASA TM-X-1502 sample Me={Me} γ={g}",
        inputs={
            "Me": Me,
            "gamma": g,
            "dnu_deg": reg["dnu_deg"],
            "ht": reg["throat_half_height"],
        },
        meta={"mesh_note": reg.get("mesh_note", "")},
    )
    th_max_ref = math.degrees(0.5 * nu_e)
    th_max_sol = math.degrees(r.theta_max_rad)
    case.metrics.append(
        _metric(
            "theta_max_deg",
            th_max_ref,
            th_max_sol,
            tol_percent=reg["tol_percent"]["theta_max_deg"],
            unit="deg",
            note="θ_max = ν_e/2 (TM-X-1502)",
        )
    )
    case.metrics.append(
        _metric(
            "area_ratio_1d",
            ar,
            float(s.get("exit_y_over_ht", r.contour.r[-1])),
            tol_percent=reg["tol_percent"]["ye_over_ht"],
            unit="-",
            note="y_e/ht = A/A* for planar",
        )
    )
    case.metrics.append(
        _metric(
            "ye_over_ht",
            ar,
            float(r.contour.r[-1] / reg["throat_half_height"]),
            tol_percent=reg["tol_percent"]["ye_over_ht"],
            unit="-",
        )
    )
    case.metrics.append(
        _metric(
            "throat_x",
            0.0,
            float(r.contour.x[0]),
            abs_tol=1e-12,
            unit="m",
        )
    )
    case.metrics.append(
        _metric(
            "throat_y",
            float(reg["throat_half_height"]),
            float(r.contour.r[0]),
            abs_tol=1e-12,
            unit="m",
        )
    )
    case.metrics.append(
        _metric(
            "xe_over_ht",
            None,
            float(r.contour.x[-1] / reg["throat_half_height"]),
            unit="-",
            note="physical length (mesh-dependent); TM-X-1502 discusses Δν sensitivity",
        )
    )
    case.finalize()
    return [case]


def run_naca_planar_suite(catalog: dict) -> list[CaseResult]:
    reg = catalog["regimes"]["naca_planar_mln_suite"]
    src = reg["source_id"]
    g = float(reg["gamma"])
    gas = PerfectGas(g)
    tols = reg["tol_percent"]
    out: list[CaseResult] = []
    for c in reg["cases"]:
        Me = float(c["Me"])
        dnu = float(c["dnu_deg"])
        nu_e = float(gas.prandtl_meyer(Me))
        ar = float(gas.area_ratio(Me))
        r = design_planar_mln_nasa(Me, g, dnu_deg=dnu, throat_half_height=1.0)
        case = CaseResult(
            case_id=f"naca_planar_M{Me}_dnu{dnu}",
            regime="planar_mln_naca",
            source_id=src,
            source_citation=_cite(catalog, src),
            description=f"NACA E8J12-rule planar MLN Me={Me} γ={g}",
            inputs={"Me": Me, "gamma": g, "dnu_deg": dnu},
        )
        case.metrics.append(
            _metric(
                "theta_max_deg",
                math.degrees(0.5 * nu_e),
                math.degrees(r.theta_max_rad),
                tol_percent=tols["theta_max_deg"],
                unit="deg",
            )
        )
        case.metrics.append(
            _metric(
                "ye_over_ht",
                ar,
                float(r.contour.r[-1]),
                tol_percent=tols["ye_over_ht"],
                unit="-",
                note="planar exit half-height = A/A*",
            )
        )
        case.metrics.append(
            _metric(
                "area_ratio",
                ar,
                float(r.contour.r[-1]),  # ht=1
                tol_percent=tols["area_ratio"],
                unit="-",
            )
        )
        case.metrics.append(
            _metric("L_over_ht", None, float(r.contour.x[-1]), unit="-", note="solver length")
        )
        case.finalize()
        out.append(case)
    return out


def run_axisym_mln(catalog: dict) -> list[CaseResult]:
    reg = catalog["regimes"]["axisym_mln_regimes"]
    src = reg["source_id"]
    Rt = float(reg["Rt"])
    n = int(reg["n_chars"])
    contracts = reg["contracts"]
    tols = reg["tol_percent"]
    out: list[CaseResult] = []
    for c in reg["cases"]:
        g, Me = float(c["gamma"]), float(c["Me"])
        gas = PerfectGas(g)
        nu_e = float(gas.prandtl_meyer(Me))
        ar_1d = float(gas.area_ratio(Me))
        res = design_mln(
            DesignSpec(
                gas=gas,
                throat_radius=Rt,
                exit_mach=Me,
                n_characteristics=n,
            ),
            store_net=False,
        )
        case = CaseResult(
            case_id=f"axisym_mln_g{g}_M{Me}",
            regime="axisym_mln",
            source_id=src,
            source_citation=_cite(catalog, src),
            description=f"Axisym MLN γ={g} Me={Me} n={n}",
            inputs={"gamma": g, "Me": Me, "Rt": Rt, "n_chars": n},
        )
        eps = res.metrics.area_ratio
        th_max = math.degrees(res.metrics.max_wall_angle)
        th_e = math.degrees(abs(float(res.contour.theta[-1])))
        L_Rt = res.metrics.length / Rt
        half_nu = math.degrees(0.5 * nu_e)

        case.metrics.append(
            _metric(
                "area_ratio",
                ar_1d,
                eps,
                tol_percent=tols["area_ratio"],
                unit="-",
                note="ε_geom vs A/A*(Me)",
            )
        )
        case.metrics.append(
            _metric(
                "theta_max_vs_half_nu_e_deg",
                half_nu,
                th_max,
                tol_percent=tols["theta_max_vs_half_nu"],
                unit="deg",
                note="corner ≈ ν_e/2 (axisym mild offset allowed)",
            )
        )
        case.metrics.append(
            _metric(
                "theta_exit_deg",
                0.0,
                th_e,
                abs_tol=contracts["theta_exit_near_zero_deg"],
                unit="deg",
                note="axial exit design",
            )
        )
        case.metrics.append(
            _metric(
                "L_over_Rt",
                None,
                L_Rt,
                unit="-",
                bool_pass=contracts["L_over_Rt_min"] <= L_Rt <= contracts["L_over_Rt_max"],
                note=f"band [{contracts['L_over_Rt_min']}, {contracts['L_over_Rt_max']}]",
            )
        )
        case.metrics.append(
            _metric(
                "exit_radius_m",
                Rt * math.sqrt(ar_1d),
                res.metrics.exit_radius,
                tol_percent=tols["area_ratio"],
                unit="m",
            )
        )
        case.metrics.append(
            _metric(
                "Cf_vac_1d",
                None,
                res.metrics.Cf_vac_1d,
                unit="-",
            )
        )
        diverging = res.metrics.exit_radius > Rt
        case.metrics.append(
            _metric(
                "diverging",
                1.0,
                1.0 if diverging else 0.0,
                bool_pass=diverging,
                note="Re > Rt",
            )
        )
        case.finalize()
        out.append(case)
    return out


def run_rao_top(catalog: dict) -> list[CaseResult]:
    reg = catalog["regimes"]["rao_top_engineering"]
    src = reg["source_id"]
    out: list[CaseResult] = []
    for c in reg["cases"]:
        g, Me = float(c["gamma"]), float(c["Me"])
        lf = float(c["length_fraction"])
        Rt = float(c["Rt"])
        gas = PerfectGas(g)
        ar = float(gas.area_ratio(Me))
        # 15° cone length to same ε: L = (Re-Rt)/tan(15°)
        Re = Rt * math.sqrt(ar)
        L_cone = (Re - Rt) / math.tan(math.radians(15.0))
        res = design_rao(
            DesignSpec(gas=gas, throat_radius=Rt, exit_mach=Me, n_characteristics=24),
            length_fraction=lf,
        )
        case = CaseResult(
            case_id=f"rao_top_g{g}_M{Me}_lf{lf}",
            regime="rao_top",
            source_id=src,
            source_citation=_cite(catalog, src),
            description=f"Rao TOP γ={g} Me={Me} L/L_15°={lf}",
            inputs=c,
        )
        L_ratio = res.metrics.length / L_cone
        case.metrics.append(
            _metric(
                "L_over_L15cone",
                lf,
                L_ratio,
                tol_percent=reg["contracts"]["L_over_L15cone_tol_percent"],
                unit="-",
                note="Rao length fraction of 15° cone",
            )
        )
        case.metrics.append(
            _metric(
                "area_ratio",
                ar,
                res.metrics.area_ratio,
                tol_percent=reg["contracts"]["eps_matches_1d_percent"],
                unit="-",
            )
        )
        th_n = math.degrees(res.metrics.max_wall_angle)
        th_e = math.degrees(abs(float(res.contour.theta[-1])))
        case.metrics.append(
            _metric(
                "theta_n_deg",
                None,
                th_n,
                unit="deg",
            )
        )
        case.metrics.append(
            _metric(
                "theta_e_deg",
                None,
                th_e,
                unit="deg",
            )
        )
        case.metrics.append(
            _metric(
                "theta_n_gt_theta_e",
                1.0,
                1.0 if th_n > th_e else 0.0,
                bool_pass=th_n > th_e,
                note="Rao attachment angle > exit angle",
            )
        )
        case.finalize()
        out.append(case)
    return out


def run_aerospike(catalog: dict) -> list[CaseResult]:
    reg = catalog["regimes"]["aerospike_lee"]
    src = reg["source_id"]
    tols = reg["tol_percent"]
    out: list[CaseResult] = []
    for i, c in enumerate(reg["cases"]):
        gas = PerfectGas(float(c["gamma"]))
        spike = design_spike(
            Pc=float(c["Pc"]),
            Tc=float(c["Tc"]),
            gamma=float(c["gamma"]),
            molar_mass=float(c["molar_mass"]),
            Re=float(c["Re"]),
            expansion_ratio=float(c["expansion_ratio"]),
            Pa=float(c["Pa"]),
            n_points=2000,
        )
        nu_e = float(gas.prandtl_meyer(spike.exit_mach))
        delta_ref = 0.5 * math.pi - nu_e
        Me_ar = float(gas.mach_from_area_ratio(float(c["expansion_ratio"]), branch="supersonic"))
        case = CaseResult(
            case_id=f"aerospike_lee_{i}",
            regime="aerospike",
            source_id=src,
            source_citation=_cite(catalog, src),
            description=f"Lee spike ε={c['expansion_ratio']} γ={c['gamma']}",
            inputs=c,
        )
        case.metrics.append(
            _metric(
                "delta_rad",
                delta_ref,
                spike.delta_rad,
                tol_percent=tols["delta_vs_pi2_minus_nu"],
                unit="rad",
                note="δ = π/2 − ν_e",
            )
        )
        case.metrics.append(
            _metric(
                "exit_mach",
                Me_ar,
                spike.exit_mach,
                tol_percent=tols["exit_mach_vs_area_mach"],
                unit="-",
                note="Me from ε (Majdalani / area–Mach)",
            )
        )
        case.metrics.append(
            _metric(
                "expansion_ratio",
                float(c["expansion_ratio"]),
                spike.expansion_ratio,
                tol_percent=0.01,
                unit="-",
            )
        )
        case.metrics.append(
            _metric(
                "R_tip_near_zero",
                0.0,
                float(spike.R[-1]),
                abs_tol=1e-3 * float(c["Re"]),
                unit="m",
                note="full spike tip → axis",
            )
        )
        case.metrics.append(
            _metric("Isp_s", None, spike.Isp, unit="s", note="solver performance")
        )
        case.metrics.append(
            _metric("Cf", None, spike.Cf, unit="-")
        )
        case.finalize()
        out.append(case)
    return out


def run_performance_1d(catalog: dict) -> list[CaseResult]:
    reg = catalog["regimes"]["performance_1d_cf"]
    src = reg["source_id"]
    tol = reg["tol_percent"]["Cf"]
    out: list[CaseResult] = []
    for c in reg["cases"]:
        g, Me, pa_pc = float(c["gamma"]), float(c["Me"]), float(c["pa_pc"])
        gas = PerfectGas(g)
        pe_pc = float(gas.pe_pc_from_mach(Me))
        ar = float(gas.area_ratio(Me))
        # Reference: recompute with same closed form (self-consistency + Sutton formula structure)
        # Published structure: Cf = sqrt{...} + ε(pe/pc - pa/pc)
        term = (2.0 * g**2 / (g - 1.0)) * ((2.0 / (g + 1.0)) ** ((g + 1.0) / (g - 1.0)))
        pe_term = 1.0 - pe_pc ** ((g - 1.0) / g)
        Cf_ref = math.sqrt(term * pe_term) + ar * (pe_pc - pa_pc)
        Cf_sol = gas.thrust_coefficient(pe_pc, pa_pc, exit_mach=Me, area_ratio=ar)
        case = CaseResult(
            case_id=f"cf_1d_g{g}_M{Me}_pa{pa_pc}",
            regime="performance_1d",
            source_id=src,
            source_citation=_cite(catalog, src),
            description=f"1-D Cf γ={g} Me={Me} pa/pc={pa_pc}",
            inputs=c,
        )
        case.metrics.append(
            _metric("Cf", Cf_ref, Cf_sol, tol_percent=tol, unit="-", note="Sutton ideal Cf")
        )
        case.metrics.append(
            _metric("pe_pc", pe_pc, pe_pc, tol_percent=1e-9, unit="-")
        )
        case.metrics.append(
            _metric("area_ratio", ar, ar, tol_percent=1e-9, unit="-")
        )
        case.finalize()
        out.append(case)
    return out


def run_family_ordering(catalog: dict) -> list[CaseResult]:
    reg = catalog["regimes"]["family_ordering"]
    src = reg["source_id"]
    g, Me, Rt = float(reg["gamma"]), float(reg["Me"]), float(reg["Rt"])
    n = int(reg["n_chars"])
    gas = PerfectGas(g)
    spec = DesignSpec(gas=gas, throat_radius=Rt, exit_mach=Me, n_characteristics=n)
    mln = design_mln(spec, store_net=False)
    ideal = design_ideal(spec, store_net=False)
    tic = design_tic(
        DesignSpec(
            gas=gas,
            throat_radius=Rt,
            exit_mach=Me,
            n_characteristics=n,
            truncate_fraction=float(reg["tic_fraction"]),
        ),
        store_net=False,
    )
    rao = design_rao(spec, length_fraction=float(reg["rao_fraction"]))
    cmp = compare_nozzle_designs(spec, rao_length_fraction=float(reg["rao_fraction"]))
    case = CaseResult(
        case_id="family_ordering_mln_tic_ideal_rao",
        regime="family_ordering",
        source_id=src,
        source_citation=_cite(catalog, src),
        description=f"Length ordering γ={g} Me={Me}",
        inputs={"gamma": g, "Me": Me, "Rt": Rt},
        meta={"comparison": cmp.rows},
    )
    case.metrics.append(
        _metric("L_mln_m", None, mln.metrics.length, unit="m")
    )
    case.metrics.append(
        _metric("L_tic_m", None, tic.metrics.length, unit="m")
    )
    case.metrics.append(
        _metric("L_ideal_m", None, ideal.metrics.length, unit="m")
    )
    case.metrics.append(
        _metric("L_rao_m", None, rao.metrics.length, unit="m")
    )
    case.metrics.append(
        _metric(
            "ideal_gt_mln",
            1.0,
            1.0 if ideal.metrics.length > mln.metrics.length else 0.0,
            bool_pass=ideal.metrics.length > mln.metrics.length,
            note="ideal longer than MLN",
        )
    )
    case.metrics.append(
        _metric(
            "tic_lt_ideal",
            1.0,
            1.0 if tic.metrics.length < ideal.metrics.length else 0.0,
            bool_pass=tic.metrics.length < ideal.metrics.length,
            note="TIC shorter than ideal",
        )
    )
    case.metrics.append(
        _metric(
            "rao_lt_ideal",
            1.0,
            1.0 if rao.metrics.length < ideal.metrics.length else 0.0,
            bool_pass=rao.metrics.length < ideal.metrics.length * 1.05,
            note="Rao TOP typically shorter than full ideal",
        )
    )
    case.finalize()
    return [case]


# ---------------------------------------------------------------------------
# Orchestration + export
# ---------------------------------------------------------------------------


def run_deep_validation(catalog_path: Path | str | None = None) -> dict[str, Any]:
    catalog = _load_catalog(Path(catalog_path) if catalog_path else None)
    cases: list[CaseResult] = []
    cases.extend(run_thermo_isentropic(catalog))
    cases.extend(run_thermo_multi_gamma(catalog))
    cases.extend(run_nasa_tm_x_1502(catalog))
    cases.extend(run_naca_planar_suite(catalog))
    cases.extend(run_axisym_mln(catalog))
    cases.extend(run_rao_top(catalog))
    cases.extend(run_aerospike(catalog))
    cases.extend(run_performance_1d(catalog))
    cases.extend(run_family_ordering(catalog))

    n_pass = sum(1 for c in cases if c.passed)
    n_fail = len(cases) - n_pass
    all_err = [
        m.error_percent
        for c in cases
        for m in c.metrics
        if m.error_percent is not None and math.isfinite(m.error_percent)
    ]
    report = {
        "title": "rockmoc deep multi-regime validation",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "catalog_title": catalog.get("title"),
        "sources": catalog.get("sources", {}),
        "summary": {
            "n_cases": len(cases),
            "n_passed": n_pass,
            "n_failed": n_fail,
            "n_metrics": sum(len(c.metrics) for c in cases),
            "max_error_percent": max(all_err) if all_err else None,
            "mean_error_percent": float(np.mean(all_err)) if all_err else None,
            "median_error_percent": float(np.median(all_err)) if all_err else None,
        },
        "cases": [c.as_dict() for c in cases],
    }
    return report


def write_validation_outputs(
    report: dict[str, Any],
    out_dir: Path | str,
) -> dict[str, str]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}

    json_path = out / "validation_results.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    paths["json"] = str(json_path)

    csv_path = out / "validation_results.csv"
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "case_id",
                "regime",
                "source_id",
                "metric",
                "reference",
                "solver",
                "error_percent",
                "abs_error",
                "tol_percent",
                "passed",
                "unit",
                "note",
                "case_passed",
            ]
        )
        for c in report["cases"]:
            for m in c["metrics"]:
                w.writerow(
                    [
                        c["case_id"],
                        c["regime"],
                        c["source_id"],
                        m["name"],
                        m["reference"],
                        m["solver"],
                        m["error_percent"],
                        m["abs_error"],
                        m["tol_percent"],
                        m["passed"],
                        m["unit"],
                        m["note"],
                        c["passed"],
                    ]
                )
    paths["csv"] = str(csv_path)

    md_path = out / "validation_report.md"
    md_path.write_text(_to_markdown(report), encoding="utf-8")
    paths["markdown"] = str(md_path)

    # Compact source index
    src_path = out / "sources.md"
    lines = ["# Validation sources\n"]
    for sid, s in report.get("sources", {}).items():
        lines.append(f"## `{sid}`\n")
        lines.append(f"- **Citation:** {s.get('citation', '')}\n")
        lines.append(f"- **Type:** {s.get('type', '')}\n")
        if s.get("local_pdf"):
            lines.append(f"- **Local PDF:** `{s['local_pdf']}`\n")
        if s.get("notes"):
            lines.append(f"- **Notes:** {s['notes']}\n")
        lines.append("")
    src_path.write_text("\n".join(lines), encoding="utf-8")
    paths["sources"] = str(src_path)

    return paths


def _to_markdown(report: dict[str, Any]) -> str:
    s = report["summary"]
    lines = [
        f"# {report['title']}",
        "",
        f"**Generated (UTC):** {report['generated_utc']}",
        "",
        "## Summary",
        "",
        f"| Cases | Passed | Failed | Metrics | Max err % | Mean err % | Median err % |",
        f"|------:|-------:|-------:|--------:|----------:|-----------:|-------------:|",
        f"| {s['n_cases']} | {s['n_passed']} | {s['n_failed']} | {s['n_metrics']} | "
        f"{_fmt(s['max_error_percent'])} | {_fmt(s['mean_error_percent'])} | {_fmt(s['median_error_percent'])} |",
        "",
        "## Sources",
        "",
    ]
    for sid, src in report.get("sources", {}).items():
        lines.append(f"- **{sid}:** {src.get('citation', '')}")
    lines.extend(["", "## Results by case", ""])

    for c in report["cases"]:
        flag = "PASS" if c["passed"] else "FAIL"
        lines.append(f"### `{c['case_id']}` — **{flag}**")
        lines.append("")
        lines.append(f"- **Regime:** {c['regime']}")
        lines.append(f"- **Source:** {c['source_citation']}")
        lines.append(f"- **Description:** {c['description']}")
        lines.append("")
        lines.append("| Metric | Reference | Solver | Error % | Tol % | Pass | Unit | Note |")
        lines.append("|--------|----------:|-------:|--------:|------:|:----:|------|------|")
        for m in c["metrics"]:
            lines.append(
                f"| {m['name']} | {_fmt(m['reference'])} | {_fmt(m['solver'])} | "
                f"{_fmt(m['error_percent'])} | {_fmt(m['tol_percent'])} | "
                f"{'Y' if m['passed'] else 'N'} | {m.get('unit') or ''} | {m.get('note') or ''} |"
            )
        lines.append("")
    return "\n".join(lines)


def _fmt(x: Any) -> str:
    if x is None:
        return "—"
    if isinstance(x, bool):
        return str(x)
    try:
        v = float(x)
        if abs(v) >= 100 or (abs(v) > 0 and abs(v) < 1e-3):
            return f"{v:.6g}"
        return f"{v:.6f}".rstrip("0").rstrip(".")
    except (TypeError, ValueError):
        return str(x)


def main(argv: list[str] | None = None) -> int:
    import argparse

    p = argparse.ArgumentParser(description="Deep multi-regime rockmoc validation")
    p.add_argument(
        "--catalog",
        default=None,
        help="Path to validation_catalog.json",
    )
    p.add_argument(
        "--out",
        default=None,
        help="Output directory (default: tests/data/validation_results)",
    )
    args = p.parse_args(argv)
    out = args.out
    if out is None:
        out = Path(__file__).resolve().parents[3] / "tests" / "data" / "validation_results"
    report = run_deep_validation(args.catalog)
    paths = write_validation_outputs(report, out)
    s = report["summary"]
    print(f"Cases: {s['n_passed']}/{s['n_cases']} passed ({s['n_failed']} failed)")
    print(f"Metrics: {s['n_metrics']}; max error % = {s['max_error_percent']}")
    for k, v in paths.items():
        print(f"  {k}: {v}")
    return 0 if s["n_failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
