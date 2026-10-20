"""Verification plots: grid convergence (Richardson / observed order)."""

from __future__ import annotations

from typing import Any

import numpy as np

from rockmoc.plot.theme import (
    COLORS,
    finalize_figure,
    maybe_save,
    maybe_show,
    require_matplotlib,
    style_axes,
    themed,
)

_LABELS = {
    "length_over_Rt": "L / $R_t$",
    "area_ratio_geometric": "ε",
    "exit_wall_angle_deg": "exit wall angle [deg]",
    "max_wall_angle_deg": "max wall angle [deg]",
}


@themed
def plot_grid_convergence(
    study: dict,
    *,
    title: str | None = None,
    show: bool = True,
    save_path: str | None = None,
    figsize: tuple[float, float] = (11, 4.5),
) -> Any:
    """Log-log error |f(N) − f_ext| vs N per metric, with slope-1 / slope-2 guides.

    ``study`` is the output of :func:`rockmoc.analysis.grid_convergence`.
    Metrics converged to round-off are skipped.
    """
    plt, _ = require_matplotlib()
    rows = study["rows"]
    n = np.array([r["n_chars"] for r in rows], dtype=float)
    stats = study["richardson"]
    metrics = [m for m, s in stats.items() if s.get("order") is not None]
    fig, axes = plt.subplots(1, max(len(metrics), 1), figsize=figsize, squeeze=False)
    fig.suptitle(title or f"Grid convergence — {study['kind']}, Me = {study['exit_mach']}",
                 color=COLORS["text"], fontweight="semibold")
    for ax, m in zip(axes[0], metrics):
        f = np.array([r[m] for r in rows])
        err = np.abs(f - stats[m]["extrapolated"])
        style_axes(ax, title=_LABELS.get(m, m), xlabel="N characteristics", ylabel="|f − f$_{ext}$|")
        ax.loglog(n[:-1], err[:-1], "o-", color=COLORS["wall"],
                  label=f"observed p = {stats[m]['order']:.2f}")
        ref = err[0] if err[0] > 0 else 1.0
        for p, ls in ((1, "--"), (2, ":")):
            ax.loglog(n, ref * (n / n[0]) ** -p, ls, color=COLORS["muted"], lw=1.0, label=f"slope −{p}")
        from matplotlib.ticker import NullFormatter, ScalarFormatter

        ax.set_xticks(n)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.legend(fontsize=8)
    finalize_figure(fig)
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig
