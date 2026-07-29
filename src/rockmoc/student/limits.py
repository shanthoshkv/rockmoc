"""Known limitations — always surface these to students and reports."""

from __future__ import annotations

KNOWN_LIMITS: list[str] = [
    "Inviscid method of characteristics — no shocks, separation, or wall heat transfer CFD.",
    "Perfect-gas or frozen tabulated γ along the nozzle — not multi-phase reacting flow.",
    "Exit performance uses 1-D / exit-station models; validate with CFD before hardware.",
    "Canceling-wall MLN is an engineering hybrid (schedule + streamline + isotropic ε match).",
    "Boundary-layer correction is an integral estimate (order-of-magnitude δ*), not full BL CFD.",
    "Aerospike uses Lee/Vernacchia isentropic annular theory (ν_e must stay below π/2).",
    "rockmoc is for preliminary design, lab reports, and contour export — not flight qualification.",
]

LIMITS_SHORT = (
    "Inviscid MOC + perfect/frozen-γ model. Preliminary design only — not a CFD substitute."
)


def limits_bullet_list(*, markdown: bool = True) -> str:
    """Format known limits for reports."""
    if markdown:
        return "\n".join(f"- {line}" for line in KNOWN_LIMITS)
    return "\n".join(f"• {line}" for line in KNOWN_LIMITS)


def limits_block() -> str:
    """Plain multi-line block for terminal reports."""
    lines = ["KNOWN LIMITS"]
    lines.append("-" * 40)
    for line in KNOWN_LIMITS:
        lines.append(f"  • {line}")
    return "\n".join(lines)
