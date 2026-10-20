"""Styling for rockmoc figures: dark (screen) and light (publication) themes.

Public plot functions are wrapped with :func:`themed`: styling is applied
inside ``matplotlib.rc_context`` so the caller's global rcParams are left
untouched, and each accepts ``theme="dark" | "light"``.
"""

from __future__ import annotations

import functools
from typing import Any

# Palette — refined dark aerospace (high contrast on OLED / projectors)
COLORS = {
    "bg": "#0b0f14",
    "panel": "#12181f",
    "panel2": "#1a222d",
    "grid": "#2a3441",
    "grid_minor": "#1e2732",
    "text": "#e8eef6",
    "muted": "#8b9aab",
    "wall": "#4da3ff",
    "wall_fill": "#1f6feb40",
    "c_plus": "#3dd68c",
    "c_minus": "#ff7b72",
    "axis": "#c4a1ff",
    "fan": "#ffa657",
    "accent": "#79c0ff",
    "warn": "#e3b341",
    "good": "#3fb950",
    "line2": "#a5d6ff",
    "line3": "#ff9b8a",
    "line4": "#d2a8ff",
    "line5": "#7ee787",
    "annotation_bg": "#0b0f14e6",
    "spine": "#3d4a5c",
}

DARK = dict(COLORS)

# Light / publication palette (print-safe, colour-blind friendly hues)
LIGHT = {
    "bg": "#ffffff",
    "panel": "#ffffff",
    "panel2": "#f6f8fa",
    "grid": "#d0d7de",
    "grid_minor": "#eaeef2",
    "text": "#1f2328",
    "muted": "#57606a",
    "wall": "#0550ae",
    "wall_fill": "#0550ae18",
    "c_plus": "#1a7f37",
    "c_minus": "#cf222e",
    "axis": "#8250df",
    "fan": "#bc4c00",
    "accent": "#0969da",
    "warn": "#9a6700",
    "good": "#1a7f37",
    "line2": "#0969da",
    "line3": "#cf222e",
    "line4": "#8250df",
    "line5": "#116329",
    "annotation_bg": "#ffffffe6",
    "spine": "#8c959f",
}

THEMES = {"dark": DARK, "light": LIGHT}


def require_matplotlib():
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required for plotting; install with: pip install rockmoc[plot]"
        ) from exc
    return plt, mpl


def apply_dark_theme(plt: Any = None, mpl: Any = None) -> None:
    """Apply the active palette (dark unless inside a light ``themed`` call) to rcParams.

    Called directly this changes the global rcParams for the session; plot
    functions call it inside ``rc_context`` via :func:`themed`.
    """
    if plt is None or mpl is None:
        plt, mpl = require_matplotlib()
    mpl.rcParams.update(
        {
            "figure.facecolor": COLORS["bg"],
            "figure.edgecolor": COLORS["bg"],
            "figure.dpi": 120,
            "savefig.dpi": 200,
            "axes.facecolor": COLORS["panel"],
            "axes.edgecolor": COLORS["spine"],
            "axes.labelcolor": COLORS["text"],
            "axes.titlecolor": COLORS["text"],
            "axes.titleweight": "semibold",
            "axes.titlepad": 12,
            "axes.labelpad": 8,
            "axes.grid": True,
            "axes.axisbelow": True,
            "axes.linewidth": 1.0,
            "grid.color": COLORS["grid"],
            "grid.linestyle": "--",
            "grid.linewidth": 0.7,
            "grid.alpha": 0.75,
            "xtick.color": COLORS["muted"],
            "ytick.color": COLORS["muted"],
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "text.color": COLORS["text"],
            "legend.facecolor": COLORS["panel2"],
            "legend.edgecolor": COLORS["spine"],
            "legend.labelcolor": COLORS["text"],
            "legend.framealpha": 0.94,
            "legend.fancybox": True,
            "legend.borderpad": 0.6,
            "legend.handlelength": 1.8,
            "savefig.facecolor": COLORS["bg"],
            "savefig.edgecolor": COLORS["bg"],
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.15,
            "font.size": 11,
            "font.family": "sans-serif",
            "axes.titlesize": 13.5,
            "axes.labelsize": 11.5,
            "figure.titlesize": 15,
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "lines.linewidth": 2.0,
            "lines.solid_capstyle": "round",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.prop_cycle": mpl.cycler(
                color=[
                    COLORS["wall"],
                    COLORS["c_plus"],
                    COLORS["c_minus"],
                    COLORS["axis"],
                    COLORS["fan"],
                    COLORS["line2"],
                    COLORS["line5"],
                ]
            ),
        }
    )


def style_axes(
    ax: Any,
    *,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    subtitle: str | None = None,
) -> None:
    """Apply panel styling, labels, and subtle polish."""
    ax.set_facecolor(COLORS["panel"])
    for spine in ax.spines.values():
        spine.set_color(COLORS["spine"])
        spine.set_linewidth(1.0)
    if title:
        ax.set_title(title, color=COLORS["text"], pad=12, fontweight="semibold", loc="left")
    if subtitle:
        ax.text(
            0.0,
            1.02,
            subtitle,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=9,
            color=COLORS["muted"],
        )
    if xlabel:
        ax.set_xlabel(xlabel, color=COLORS["text"])
    if ylabel:
        ax.set_ylabel(ylabel, color=COLORS["text"])
    ax.tick_params(colors=COLORS["muted"], which="both")
    ax.grid(True, color=COLORS["grid"], linestyle="--", alpha=0.75, which="major")
    ax.set_axisbelow(True)


def annotation_box(
    ax: Any,
    text: str,
    *,
    x: float = 0.02,
    y: float = 0.98,
    fontsize: float = 9,
) -> None:
    """Semi-transparent metrics / note box in axes coordinates."""
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=fontsize,
        color=COLORS["text"],
        family="monospace",
        linespacing=1.35,
        bbox=dict(
            boxstyle="round,pad=0.45",
            facecolor=COLORS["annotation_bg"],
            edgecolor=COLORS["spine"],
            alpha=0.95,
            linewidth=0.9,
        ),
        zorder=20,
    )


def style_colorbar(cbar: Any, label: str = "") -> None:
    cbar.set_label(label, color=COLORS["text"])
    cbar.ax.yaxis.set_tick_params(color=COLORS["muted"])
    cbar.outline.set_edgecolor(COLORS["spine"])  # type: ignore[attr-defined]
    plt, _ = require_matplotlib()
    plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color=COLORS["muted"])


def watermark(fig: Any, text: str = "rockmoc") -> None:
    fig.text(
        0.99,
        0.01,
        text,
        ha="right",
        va="bottom",
        fontsize=8,
        color=COLORS["muted"],
        alpha=0.55,
    )


def finalize_figure(fig: Any) -> None:
    """Tight layout + watermark."""
    try:
        fig.tight_layout()
    except Exception:
        pass
    watermark(fig)


def maybe_save(fig: Any, path: str | None, *, dpi: int = 200) -> str | None:
    """Save ``fig``. No suffix → SVG (vector, LaTeX/Inkscape friendly).

    ``.svg`` / ``.pdf`` / ``.eps`` are written as vectors (dpi only affects
    embedded rasters such as filled contours); other suffixes as rasters.
    """
    if path is None:
        return None
    from pathlib import Path

    from rockmoc._version import __version__

    p = Path(path)
    if not p.suffix:
        p = p.with_suffix(".svg")
    p.parent.mkdir(parents=True, exist_ok=True)
    meta = {"Creator": f"rockmoc {__version__}"} if p.suffix.lower() in (".svg", ".pdf") else None
    fig.savefig(p, dpi=dpi, bbox_inches="tight", facecolor=COLORS["bg"], edgecolor=COLORS["bg"],
                metadata=meta)
    return str(p)


def themed(fn: Any) -> Any:
    """Decorator: run a plot function with ``theme`` palette inside ``rc_context``."""

    @functools.wraps(fn)
    def wrapper(*args: Any, theme: str = "dark", **kwargs: Any) -> Any:
        if theme not in THEMES:
            raise ValueError(f"theme must be one of {sorted(THEMES)}")
        plt, mpl = require_matplotlib()
        saved = dict(COLORS)
        COLORS.update(THEMES[theme])
        try:
            with mpl.rc_context():
                apply_dark_theme(plt, mpl)
                return fn(*args, **kwargs)
        finally:
            COLORS.clear()
            COLORS.update(saved)

    return wrapper


def maybe_show(fig: Any, show: bool) -> None:
    if show:
        plt, _ = require_matplotlib()
        plt.show()
    else:
        plt, _ = require_matplotlib()
        plt.close(fig)
