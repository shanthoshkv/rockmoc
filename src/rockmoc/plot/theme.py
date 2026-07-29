"""Dark-theme styling for rockmoc figures (publication-oriented)."""

from __future__ import annotations

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
    if path is None:
        return None
    from pathlib import Path

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=dpi, bbox_inches="tight", facecolor=COLORS["bg"], edgecolor=COLORS["bg"])
    return str(p)


def maybe_show(fig: Any, show: bool) -> None:
    if show:
        plt, _ = require_matplotlib()
        plt.show()
    else:
        plt, _ = require_matplotlib()
        plt.close(fig)
