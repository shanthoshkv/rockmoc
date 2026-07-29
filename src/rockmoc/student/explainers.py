"""Short “why this shape?” explainers for each nozzle family."""

from __future__ import annotations

EXPLAINERS: dict[str, dict[str, object]] = {
    "mln": {
        "name": "Minimum-Length Nozzle (MLN)",
        "one_liner": "Shortest contoured nozzle for a given design Me (sharp throat + canceling wall).",
        "bullets": [
            "Sharp throat corner starts a Prandtl–Meyer expansion fan.",
            "Corner angle θ_max ≈ ν_e/2 (planar rule; mild axisym correction).",
            "Canceling wall turns the flow until θ_e ≈ 0 at design Me.",
            "Best when length/mass is critical (upper stages, high ε vacuum).",
            "Exit is design-uniform only approximately (MOC engineering model).",
        ],
        "use_when": "You want the shortest high-ε bell for a target Mach.",
        "compare_to": "Shorter than ideal and typical TIC; often shorter than Rao TOP.",
    },
    "ideal": {
        "name": "Ideal contoured bell",
        "one_liner": "Longer classical contoured nozzle: throat arc + canceling region toward uniform exit.",
        "bullets": [
            "Circular-arc throat (Rc/Rt) then elongated expansion/canceling region.",
            "Historically the reference for “perfect” inviscid exit uniformity.",
            "Longer and heavier than MLN or truncated designs.",
            "Useful as an upper-bound length reference.",
        ],
        "use_when": "You need a full-length contoured reference or textbook ideal shape.",
        "compare_to": "Longer than MLN, TIC, and typical Rao TOP at the same Me.",
    },
    "tic": {
        "name": "Truncated Ideal Contour (TIC)",
        "one_liner": "Cut an ideal contour to a fraction of its length to save mass.",
        "bullets": [
            "Start from ideal MOC/contour family, truncate axially (e.g. 80%).",
            "Exit angle θ_e is usually nonzero → divergence loss vs full ideal.",
            "Classic length–performance trade used on many flight engines.",
            "truncate_fraction is of ideal length (0.05 < f ≤ 1).",
        ],
        "use_when": "You like the ideal family but cannot afford full length.",
        "compare_to": "Shorter than ideal; longer or shorter than MLN depending on Me and f.",
    },
    "rao": {
        "name": "Rao TOP (thrust-optimized parabola)",
        "one_liner": "Industry engineering contour: high thrust for a fixed length (often ~80% of 15° cone).",
        "bullets": [
            "Parabolic / TOP correlations approximate Rao’s optimum-thrust idea.",
            "length_fraction is relative to a 15° conical nozzle of the same ε.",
            "Attachment angle θ_n > exit angle θ_e by design.",
            "Default choice for many liquid-engine preliminary designs.",
        ],
        "use_when": "You want a compact high-thrust bell used in industry practice.",
        "compare_to": "Usually shorter than full ideal; not identical to pure MOC MLN.",
    },
    "conical": {
        "name": "Conical nozzle",
        "one_liner": "Simplest diverging cone (default half-angle 15°) — textbook baseline.",
        "bullets": [
            "Straight wall at constant half-angle θ_c.",
            "Length L = (R_e − R_t) / tan(θ_c) with R_e from 1-D ε(Me).",
            "Divergence efficiency η_div ≈ (1 + cos θ_c)/2 (not axial exit).",
            "Use as a baseline before claiming MOC/TOP gains.",
        ],
        "use_when": "Homework baseline or “is MOC worth it?” comparisons.",
        "compare_to": "Usually longer than Rao TOP at 15° for the same ε; large θ_e.",
    },
    "planar_mln": {
        "name": "Planar MLN (NASA TM-X-1502)",
        "one_liner": "2D sharp-throat minimum-length nozzle with uniform exit (paper-faithful).",
        "bullets": [
            "θ_max = ν_e/2 exactly (planar theory).",
            "Exit half-height y_e/y_t = A/A*(Me).",
            "Mesh controlled by Δν in degrees (finer → longer runtime).",
            "Great for learning MOC and matching published NASA samples.",
        ],
        "use_when": "Course labs on planar MOC or TM-X-1502 reproduction.",
        "compare_to": "Not axisymmetric — do not compare L/Rt directly to bell nozzles.",
    },
    "aerospike": {
        "name": "Annular aerospike (Lee / Vernacchia)",
        "one_liner": "Altitude-compensating plug nozzle from isentropic spike theory.",
        "bullets": [
            "Spike contour from Prandtl–Meyer turning; δ = π/2 − ν_e.",
            "Open wake / ambient pressure affects performance with altitude.",
            "Truncation + base pressure models for practical plugs.",
            "Requires ν_e < π/2 (very high ε + high γ can fail).",
        ],
        "use_when": "Exploring altitude compensation or plug geometry.",
        "compare_to": "Different architecture than closed bells — compare Isp vs altitude, not L/Rt alone.",
    },
    "production": {
        "name": "Production design package",
        "one_liner": "Fine MOC + optional multi-γ + BL δ* + field Cf + altitude sweep.",
        "bullets": [
            "Wraps MLN/ideal/TIC/Rao with accuracy extras.",
            "BL correction grows the wall by displacement thickness.",
            "Field Cf integrates exit-plane momentum/pressure.",
            "Best “submit this contour” path for projects.",
        ],
        "use_when": "Final lab package or export for CFD.",
        "compare_to": "Same family as underlying kind, with viscous/performance post-processing.",
    },
}


def explainer(kind: str) -> str:
    """Human-readable explainer for a nozzle kind key."""
    key = kind.lower().strip()
    if key not in EXPLAINERS:
        known = ", ".join(sorted(EXPLAINERS))
        return f"No explainer for {kind!r}. Known: {known}."
    e = EXPLAINERS[key]
    lines = [
        f"{e['name']}",
        f"  {e['one_liner']}",
        "  Why this shape:",
    ]
    for b in e["bullets"]:  # type: ignore[union-attr]
        lines.append(f"    • {b}")
    lines.append(f"  Use when: {e['use_when']}")
    lines.append(f"  Compare to: {e['compare_to']}")
    return "\n".join(lines)


def explainer_markdown(kind: str) -> str:
    """Markdown version of the explainer."""
    key = kind.lower().strip()
    if key not in EXPLAINERS:
        return f"_No explainer for `{kind}`._"
    e = EXPLAINERS[key]
    lines = [
        f"### {e['name']}",
        "",
        e["one_liner"],  # type: ignore[list-item]
        "",
        "**Why this shape**",
        "",
    ]
    for b in e["bullets"]:  # type: ignore[union-attr]
        lines.append(f"- {b}")
    lines.extend(
        [
            "",
            f"**Use when:** {e['use_when']}",
            "",
            f"**Compare to:** {e['compare_to']}",
            "",
        ]
    )
    return "\n".join(lines)


def all_explainers_markdown() -> str:
    """Full mini-explainer sheet for docs."""
    parts = ["# Nozzle type mini-explainers", ""]
    for key in ("mln", "ideal", "tic", "rao", "conical", "planar_mln", "aerospike", "production"):
        parts.append(explainer_markdown(key))
    return "\n".join(parts)
