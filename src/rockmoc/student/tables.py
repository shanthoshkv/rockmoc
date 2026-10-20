"""Side-by-side comparison tables (Markdown / LaTeX / CSV / plain)."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Sequence

_DEFAULT_COLUMNS = (
    ("kind", "kind"),
    ("L_over_Rt", "L/Rt"),
    ("area_ratio", "ε"),
    ("Cf_vac", "Cf,vac"),
    ("Isp_vac_s", "Isp,vac [s]"),
    ("eta_n", "η_n"),
    ("theta_max_deg", "θ_max [deg]"),
    ("theta_exit_deg", "θ_e [deg]"),
)


def _rows_from_comparison(comparison: Any) -> list[dict[str, Any]]:
    if hasattr(comparison, "rows"):
        return list(comparison.rows)
    if isinstance(comparison, list):
        return comparison
    raise TypeError("Expected DesignComparison or list[dict] rows")


def comparison_to_markdown(
    comparison: Any,
    *,
    columns: Sequence[tuple[str, str]] | None = None,
    title: str | None = "Nozzle family comparison",
) -> str:
    """GitHub-flavored Markdown table."""
    rows = _rows_from_comparison(comparison)
    cols = list(columns) if columns is not None else list(_DEFAULT_COLUMNS)
    # drop missing keys
    if rows:
        cols = [(k, h) for k, h in cols if any(k in r for r in rows)]
    headers = [h for _, h in cols]
    lines: list[str] = []
    if title:
        lines.extend([f"### {title}", ""])
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for r in rows:
        cells = []
        for key, _h in cols:
            v = r.get(key, "")
            if isinstance(v, float):
                cells.append(f"{v:.4g}")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")
    return "\n".join(lines)


def comparison_to_latex(
    comparison: Any,
    *,
    columns: Sequence[tuple[str, str]] | None = None,
    caption: str = "Nozzle family comparison",
    label: str = "tab:rockmoc_compare",
) -> str:
    """Simple LaTeX tabular (booktabs-friendly)."""
    rows = _rows_from_comparison(comparison)
    cols = list(columns) if columns is not None else list(_DEFAULT_COLUMNS)
    if rows:
        cols = [(k, h) for k, h in cols if any(k in r for r in rows)]
    n = len(cols)
    align = "l" + "r" * (n - 1)
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        rf"\caption{{{caption}}}",
        rf"\label{{{label}}}",
        rf"\begin{{tabular}}{{{align}}}",
        r"\toprule",
        " & ".join(_latex_escape(h) for _, h in cols) + r" \\",
        r"\midrule",
    ]
    for r in rows:
        cells = []
        for key, _h in cols:
            v = r.get(key, "")
            if isinstance(v, float):
                cells.append(f"{v:.4g}")
            else:
                cells.append(_latex_escape(str(v)))
        lines.append(" & ".join(cells) + r" \\")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table}",
            "",
        ]
    )
    return "\n".join(lines)


def comparison_to_csv(
    comparison: Any,
    path: str | Path | None = None,
    *,
    columns: Sequence[tuple[str, str]] | None = None,
) -> str:
    """Write or return CSV text of comparison rows."""
    rows = _rows_from_comparison(comparison)
    cols = list(columns) if columns is not None else list(_DEFAULT_COLUMNS)
    if rows:
        # include all keys present
        keys = []
        seen = set()
        for k, _h in cols:
            if k not in seen and any(k in r for r in rows):
                keys.append(k)
                seen.add(k)
        for r in rows:
            for k in r:
                if k not in seen:
                    keys.append(k)
                    seen.add(k)
    else:
        keys = [k for k, _ in cols]

    import io

    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in keys})
    text = buf.getvalue()
    if path is not None:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    return text


def comparison_to_string(comparison: Any) -> str:
    """Plain fixed-width table (uses DesignComparison.summary_table if present)."""
    if hasattr(comparison, "summary_table"):
        return comparison.summary_table()
    return comparison_to_markdown(comparison, title=None)


def _latex_escape(s: str) -> str:
    return (
        s.replace("\\", r"\textbackslash{}")
        .replace("&", r"\&")
        .replace("%", r"\%")
        .replace("_", r"\_")
        .replace("#", r"\#")
    )
