"""Unit-aware design reports and lab-report checklists for students."""

from __future__ import annotations

import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from rockmoc.student.explainers import explainer, explainer_markdown
from rockmoc.student.limits import KNOWN_LIMITS, LIMITS_SHORT, limits_bullet_list

LengthUnit = Literal["m", "mm", "cm", "in"]
AngleUnit = Literal["rad", "deg"]


def _scale_length(meters: float, unit: LengthUnit) -> tuple[float, str]:
    if unit == "m":
        return meters, "m"
    if unit == "mm":
        return meters * 1000.0, "mm"
    if unit == "cm":
        return meters * 100.0, "cm"
    if unit == "in":
        return meters / 0.0254, "in"
    raise ValueError(f"Unknown length unit {unit!r}")


def _fmt_angle(rad: float, unit: AngleUnit) -> tuple[float, str]:
    if unit == "rad":
        return rad, "rad"
    return math.degrees(rad), "deg"


def format_design_report(
    result: Any,
    *,
    length_unit: LengthUnit = "mm",
    angle_unit: AngleUnit = "deg",
    include_limits: bool = True,
    include_explainer: bool = True,
    title: str | None = None,
) -> str:
    """Human-readable lab-style report for a DesignResult (or compatible)."""
    kind = getattr(result, "kind", "design")
    spec = result.spec
    m = result.metrics
    Rt = spec.throat_radius
    Me = float(result.exit_mach_design)
    g = float(spec.gas.gamma)

    Rt_v, lu = _scale_length(Rt, length_unit)
    Re_v, _ = _scale_length(m.exit_radius, length_unit)
    L_v, _ = _scale_length(m.length, length_unit)
    th_max, au = _fmt_angle(m.max_wall_angle, angle_unit)
    th_e_rad = float(result.contour.theta[-1]) if result.contour.theta is not None else 0.0
    th_e, _ = _fmt_angle(abs(th_e_rad), angle_unit)

    lines = [
        title or "=== rockmoc design report ===",
        f"Type: {str(kind).upper()}",
        f"γ = {g:.4g}    R = {spec.gas.R:.4g} J/(kg·K)",
        f"Rt = {Rt_v:.4g} {lu}    Me = {Me:.4g}",
        f"ε = {m.area_ratio:.4g}    Re = {Re_v:.4g} {lu}",
        f"L = {L_v:.4g} {lu}    L/Rt = {m.length / Rt:.4g}",
        f"θ_max = {th_max:.3g} {au}    θ_e = {th_e:.3g} {au}",
        f"Cf,vac (1-D) = {m.Cf_vac_1d:.4g}    pe/pc = {m.pe_pc_1d:.4g}",
        f"n_characteristics = {spec.n_characteristics}",
    ]
    if include_explainer:
        key = str(kind)
        if key in ("mln", "ideal", "tic", "rao", "conical"):
            lines.append("")
            lines.append(explainer(key))
    if include_limits:
        lines.append("")
        lines.append(f"Limits: {LIMITS_SHORT}")
    lines.append("")
    lines.append("Next: export_csv → CAD/CFD; double n once to check mesh sensitivity.")
    return "\n".join(lines)


def design_report_dict(
    result: Any,
    *,
    length_unit: LengthUnit = "mm",
    angle_unit: AngleUnit = "deg",
) -> dict[str, Any]:
    """Structured unit-aware summary for JSON / notebooks."""
    m = result.metrics
    Rt = result.spec.throat_radius
    Rt_v, lu = _scale_length(Rt, length_unit)
    Re_v, _ = _scale_length(m.exit_radius, length_unit)
    L_v, _ = _scale_length(m.length, length_unit)
    th_max, au = _fmt_angle(m.max_wall_angle, angle_unit)
    th_e_rad = float(result.contour.theta[-1]) if result.contour.theta is not None else 0.0
    th_e, _ = _fmt_angle(abs(th_e_rad), angle_unit)
    return {
        "kind": result.kind,
        "gamma": result.spec.gas.gamma,
        "R_J_per_kgK": result.spec.gas.R,
        "throat_radius": Rt_v,
        "exit_radius": Re_v,
        "length": L_v,
        "length_unit": lu,
        "exit_mach": result.exit_mach_design,
        "area_ratio": m.area_ratio,
        "L_over_Rt": m.length / Rt,
        "theta_max": th_max,
        "theta_exit": th_e,
        "angle_unit": au,
        "Cf_vac_1d": m.Cf_vac_1d,
        "pe_pc_1d": m.pe_pc_1d,
        "n_characteristics": result.spec.n_characteristics,
        "limitations": LIMITS_SHORT,
        "known_limits": list(KNOWN_LIMITS),
    }


def generate_lab_report(
    result: Any,
    path: str | Path,
    *,
    student_name: str = "",
    course: str = "",
    length_unit: LengthUnit = "mm",
    angle_unit: AngleUnit = "deg",
    include_checklist: bool = True,
    extra_notes: str = "",
) -> str:
    """Write a Markdown lab-report skeleton with filled metrics + checklist.

    Returns the path written.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    kind = str(result.kind)
    d = design_report_dict(result, length_unit=length_unit, angle_unit=angle_unit)
    lu = d["length_unit"]
    au = d["angle_unit"]

    lines = [
        f"# Nozzle design lab report — {kind.upper()}",
        "",
        f"**Generated (UTC):** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M')}",
        f"**Student:** {student_name or '_your name_'}",
        f"**Course:** {course or '_course / section_'}",
        f"**Tool:** rockmoc",
        "",
        "## 1. Objective",
        "",
        f"Design a **{kind}** rocket nozzle contour for Me = {d['exit_mach']:.4g}, "
        f"γ = {d['gamma']:.4g}, and report geometry and 1-D performance.",
        "",
        "## 2. Inputs",
        "",
        "| Quantity | Value |",
        "|----------|------:|",
        f"| Kind | `{kind}` |",
        f"| γ | {d['gamma']:.6g} |",
        f"| R | {d['R_J_per_kgK']:.6g} J/(kg·K) |",
        f"| Throat radius Rt | {d['throat_radius']:.6g} {lu} |",
        f"| Exit Mach Me | {d['exit_mach']:.6g} |",
        f"| n_characteristics | {d['n_characteristics']} |",
        "",
        "## 3. Theory notes (fill in your words)",
        "",
        explainer_markdown(kind if kind in ("mln", "ideal", "tic", "rao", "conical") else "mln"),
        "",
        "Key relations used by rockmoc (perfect gas):",
        "",
        r"- Area–Mach: $A/A^* = f(M,\gamma)$",
        r"- Prandtl–Meyer: $\nu = \nu(M,\gamma)$",
        r"- Planar MLN corner: $\theta_{\max} = \nu_e/2$",
        r"- Ideal vacuum $C_F$ from isentropic exit state + geometric $\varepsilon$",
        "",
        "## 4. Results",
        "",
        "| Metric | Value |",
        "|--------|------:|",
        f"| Area ratio ε | {d['area_ratio']:.6g} |",
        f"| Exit radius Re | {d['exit_radius']:.6g} {lu} |",
        f"| Length L | {d['length']:.6g} {lu} |",
        f"| L / Rt | {d['L_over_Rt']:.6g} |",
        f"| θ_max | {d['theta_max']:.6g} {au} |",
        f"| θ_e | {d['theta_exit']:.6g} {au} |",
        f"| Cf,vac (1-D) | {d['Cf_vac_1d']:.6g} |",
        f"| pe / pc | {d['pe_pc_1d']:.6g} |",
        "",
        "```",
        format_design_report(
            result,
            length_unit=length_unit,
            angle_unit=angle_unit,
            include_limits=False,
            include_explainer=False,
        ),
        "```",
        "",
        "## 5. Figures",
        "",
        "_Attach plots from `result.plot(save_dir=...)` or notebook cells:_",
        "",
        "- [ ] Wall contour (mirrored)",
        "- [ ] MOC field (if `store_net=True`)",
        "- [ ] Optional: family comparison or altitude sweep",
        "",
        "## 6. Known limitations",
        "",
        limits_bullet_list(markdown=True),
        "",
    ]
    if include_checklist:
        lines.extend(
            [
                "## 7. Student checklist",
                "",
                "- [ ] Stated Me, γ, Rt (with units) in the objective",
                "- [ ] Explained why this nozzle type was chosen",
                "- [ ] Reported ε, L, L/Rt, θ_max, θ_e",
                "- [ ] Included at least one figure of the wall",
                "- [ ] Discussed mesh: what `n_characteristics` was used; optional double-n check",
                "- [ ] Compared to a **conical 15°** baseline or another family (table)",
                "- [ ] Listed limitations (inviscid MOC, frozen γ, not CFD)",
                "- [ ] Exported CSV path noted for CAD/CFD handoff",
                "- [ ] Cited rockmoc version and any textbooks/NASA reports used",
                "",
            ]
        )
    if extra_notes:
        lines.extend(["## 8. Additional notes", "", extra_notes, ""])
    lines.extend(
        [
            "---",
            "",
            "*Auto-generated skeleton by `rockmoc.student.generate_lab_report` — complete the narrative yourself.*",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")
    return str(path)
