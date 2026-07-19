"""Dark-theme styling for rockmoc figures."""

from __future__ import annotations

from typing import Any

# Palette (dark aerospace)
COLORS = {
    "bg": "#0d1117",
    "panel": "#161b22",
    "grid": "#30363d",
    "text": "#e6edf3",
    "muted": "#8b949e",
    "wall": "#58a6ff",
    "wall_fill": "#1f6feb33",
    "c_plus": "#3fb950",
    "c_minus": "#f78166",
    "axis": "#d2a8ff",
    "fan": "#ffa657",
    "accent": "#79c0ff",
    "warn": "#f2cc60",
    "good": "#3fb950",
    "line2": "#a5d6ff",
    "line3": "#ff7b72",
    "line4": "#d2a8ff",
    "line5": "#7ee787",
}


def require_matplotlib():
    try:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required for plotting; install with: pip install rockmoc[plot]"
        ) from exc
    return plt, mpl


def apply_dark_theme(plt: Any = None, mpl: Any = None) -> None:
    """Apply a consistent dark theme globally for the session."""
    if plt is None or mpl is None:
        plt, mpl = require_matplotlib()
    mpl.rcParams.update(
        {
            "figure.facecolor": COLORS["bg"],
            "axes.facecolor": COLORS["panel"],
            "axes.edgecolor": COLORS["grid"],
            "axes.labelcolor": COLORS["text"],
            "axes.titlecolor": COLORS["text"],
            "axes.grid": True,
            "grid.color": COLORS["grid"],
            "grid.linestyle": "--",
            "grid.linewidth": 0.6,
            "grid.alpha": 0.85,
            "xtick.color": COLORS["muted"],
            "ytick.color": COLORS["muted"],
            "text.color": COLORS["text"],
            "legend.facecolor": COLORS["panel"],
            "legend.edgecolor": COLORS["grid"],
            "legend.labelcolor": COLORS["text"],
            "legend.framealpha": 0.92,
            "savefig.facecolor": COLORS["bg"],
            "savefig.edgecolor": COLORS["bg"],
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "figure.titlesize": 14,
            "lines.linewidth": 1.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def style_axes(ax: Any, *, title: str | None = None, xlabel: str | None = None, ylabel: str | None = None) -> None:
    ax.set_facecolor(COLORS["panel"])
    for spine in ax.spines.values():
        spine.set_color(COLORS["grid"])
    if title:
        ax.set_title(title, color=COLORS["text"], pad=10, fontweight="semibold")
    if xlabel:
        ax.set_xlabel(xlabel, color=COLORS["text"])
    if ylabel:
        ax.set_ylabel(ylabel, color=COLORS["text"])
    ax.tick_params(colors=COLORS["muted"])
    ax.grid(True, color=COLORS["grid"], linestyle="--", alpha=0.8)


def maybe_save(fig: Any, path: str | None, *, dpi: int = 160) -> str | None:
    if path is None:
        return None
    from pathlib import Path

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=dpi, bbox_inches="tight", facecolor=COLORS["bg"])
    return str(p)


def maybe_show(fig: Any, show: bool) -> None:
    if show:
        plt, _ = require_matplotlib()
        plt.show()
    else:
        plt, _ = require_matplotlib()
        plt.close(fig)
