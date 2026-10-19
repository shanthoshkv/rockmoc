"""MOC characteristic-field plots (C+, C−, wall, axis) — dark theme.

Draws the **full** characteristic net when ``net.interior`` rows span the
nozzle (axis → wall). Fan rays, C− (within-row), C+ (across-row), wall
nodes, and the nozzle contour are layered for textbook-style MOC figures.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from rockmoc.geometry.contour import Contour
from rockmoc.moc.types import CharNet, FlowPoint
from rockmoc.plot.theme import (
    COLORS,
    annotation_box,
    apply_dark_theme,
    finalize_figure,
    maybe_save,
    maybe_show,
    require_matplotlib,
    style_axes,
)


def _structured_rows(net: CharNet) -> list[list[FlowPoint]]:
    """Prefer net.interior rows sorted by radius (axis → wall)."""
    if not net.interior:
        return []
    rows = []
    for row in net.interior:
        pts = [p for p in row if p is not None]
        if len(pts) >= 2:
            pts = sorted(pts, key=lambda p: (p.r, p.x))
            rows.append(pts)
    return rows


def _plot_c_families(
    ax: Any,
    net: CharNet,
    *,
    mirror: bool = True,
    alpha: float = 0.85,
    max_lines: int = 80,
) -> None:
    """Draw C− (within row) and C+ (across rows) from structured kernel."""
    rows = _structured_rows(net)

    if rows:
        # Subsample dense meshes for readability
        step = max(1, len(rows) // max_lines)
        # C− family: within each row (axis → wall) — left-running characteristics
        labeled_m = False
        for i, row in enumerate(rows):
            if i % step != 0 and i != len(rows) - 1:
                continue
            xs = [p.x for p in row]
            rs = [p.r for p in row]
            lab = "C− family" if not labeled_m else None
            ax.plot(
                xs, rs,
                color=COLORS["c_minus"],
                lw=1.05,
                alpha=alpha,
                label=lab,
                zorder=3,
            )
            if mirror:
                ax.plot(
                    xs, [-r for r in rs],
                    color=COLORS["c_minus"],
                    lw=0.75,
                    alpha=alpha * 0.35,
                    zorder=2,
                )
            labeled_m = True

        # C+ family: connect equal index across rows (right-running)
        max_j = max(len(r) for r in rows)
        step_j = max(1, max_j // max_lines)
        labeled_p = False
        for j in range(0, max_j, step_j):
            xs, rs = [], []
            for row in rows:
                if j < len(row):
                    xs.append(row[j].x)
                    rs.append(row[j].r)
            if len(xs) < 2:
                continue
            lab = "C+ family" if not labeled_p else None
            ax.plot(
                xs, rs,
                color=COLORS["c_plus"],
                lw=1.05,
                alpha=alpha,
                label=lab,
                zorder=3,
            )
            if mirror:
                ax.plot(
                    xs, [-r for r in rs],
                    color=COLORS["c_plus"],
                    lw=0.75,
                    alpha=alpha * 0.35,
                    zorder=2,
                )
            labeled_p = True
        return

    # Fallback: scatter + Riemann bins
    if not net.points:
        return
    pts = net.points
    ax.plot(
        [p.x for p in pts],
        [p.r for p in pts],
        ".",
        ms=2.5,
        color=COLORS["accent"],
        alpha=0.4,
        label="MOC nodes",
        zorder=3,
    )
    if mirror:
        ax.plot(
            [p.x for p in pts],
            [-p.r for p in pts],
            ".",
            ms=2.0,
            color=COLORS["accent"],
            alpha=0.15,
            zorder=2,
        )


def _plot_fan_rays(
    ax: Any,
    net: CharNet,
    *,
    Rt: float,
    mirror: bool = True,
) -> None:
    """Draw expansion-fan rays from the throat corner to first-wave nodes."""
    corner_x, corner_r = 0.0, Rt
    # Prefer first interior row (wave 0) points excluding pure axis
    targets: list[FlowPoint] = []
    if net.interior:
        row0 = net.interior[0]
        targets = [p for p in row0 if p.r > 1e-9]
    elif net.fan:
        targets = list(net.fan)

    if not targets:
        return
    labeled = False
    # Subsample
    step = max(1, len(targets) // 24)
    for i, p in enumerate(targets):
        if i % step != 0 and i != len(targets) - 1:
            continue
        lab = "Expansion fan" if not labeled else None
        ax.plot(
            [corner_x, p.x],
            [corner_r, p.r],
            color=COLORS["fan"],
            lw=0.9,
            alpha=0.55,
            label=lab,
            zorder=4,
        )
        if mirror:
            ax.plot(
                [corner_x, p.x],
                [-corner_r, -p.r],
                color=COLORS["fan"],
                lw=0.7,
                alpha=0.25,
                zorder=3,
            )
        labeled = True
    # Fan markers
    if targets:
        ax.plot(
            [p.x for p in targets[::step]],
            [p.r for p in targets[::step]],
            "o",
            ms=3.0,
            color=COLORS["fan"],
            alpha=0.85,
            zorder=6,
        )


def _plot_wall_nodes(ax: Any, wall: list[FlowPoint], *, mirror: bool = True) -> None:
    if not wall:
        return
    xs = [p.x for p in wall]
    rs = [p.r for p in wall]
    ax.plot(
        xs, rs,
        "o",
        ms=3.2,
        color=COLORS["wall"],
        alpha=0.9,
        label="Wall nodes",
        zorder=6,
        markeredgecolor=COLORS.get("bg", "#0b0f14"),
        markeredgewidth=0.4,
    )
    if mirror:
        ax.plot(xs, [-r for r in rs], "o", ms=2.2, color=COLORS["wall"], alpha=0.35, zorder=5)


def plot_moc_field(
    contour: Contour,
    net: CharNet | None = None,
    *,
    title: str = "Method of Characteristics — field & wall",
    mirror: bool = True,
    show: bool = True,
    save_path: str | None = None,
    figsize: tuple[float, float] = (11, 6.5),
    annotate: dict | None = None,
    show_canceling: bool = True,
) -> Any:
    """Dark-themed MOC plot: wall, axis, full C+/C− mesh, fan rays.

    Parameters
    ----------
    show_canceling:
        Kept for API compatibility. Full nets already include canceling-region
        characteristics in ``net.interior``; no synthetic stubs are drawn.
    """
    _ = show_canceling  # full net supersedes stub canceling rays
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)

    fig, ax = plt.subplots(figsize=figsize)
    style_axes(ax, title=title, xlabel="Axial station  x  [m]", ylabel="Radius  r  [m]")

    Rt = float(contour.r[0]) if len(contour.r) else 0.05

    # Characteristics under the wall
    if net is not None:
        _plot_c_families(ax, net, mirror=mirror)
        _plot_fan_rays(ax, net, Rt=Rt, mirror=mirror)

        if net.axis:
            ax.plot(
                [p.x for p in net.axis],
                [p.r for p in net.axis],
                color=COLORS["axis"],
                lw=1.6,
                ls=":",
                marker=".",
                ms=4,
                label="Axis nodes",
                zorder=4,
            )

        if net.wall:
            _plot_wall_nodes(ax, net.wall, mirror=mirror)

    # Nozzle wall (upper + mirror)
    ax.plot(contour.x, contour.r, color=COLORS["wall"], lw=2.9, label="Nozzle wall", zorder=7)
    if mirror:
        ax.plot(contour.x, -contour.r, color=COLORS["wall"], lw=2.9, zorder=7)
        ax.fill_between(
            contour.x,
            contour.r,
            -contour.r,
            color=COLORS["wall_fill"],
            alpha=0.20,
            zorder=1,
            label="Flow passage",
        )

    # Centerline
    x0, x1 = float(contour.x[0]), float(contour.x[-1])
    ax.axhline(0.0, color=COLORS["muted"], lw=1.0, ls="-", alpha=0.7, label="Centerline", zorder=1)

    # Throat / exit markers
    ax.plot([contour.x[0]], [contour.r[0]], "s", color=COLORS["warn"], ms=9, zorder=8, label="Throat")
    ax.plot([contour.x[-1]], [contour.r[-1]], "D", color=COLORS["good"], ms=9, zorder=8, label="Exit")
    if mirror:
        ax.plot([contour.x[0]], [-contour.r[0]], "s", color=COLORS["warn"], ms=7, zorder=8)
        ax.plot([contour.x[-1]], [-contour.r[-1]], "D", color=COLORS["good"], ms=7, zorder=8)

    ax.set_aspect("equal", adjustable="datalim")
    rmax = float(np.max(np.abs(contour.r))) * 1.18 + 1e-6
    xmin = x0
    xmax = x1
    if net is not None and net.points:
        xmin = min(xmin, min(p.x for p in net.points))
        xmax = max(xmax, max(p.x for p in net.points))
    pad = 0.05 * (xmax - xmin + 1e-9)
    ax.set_ylim(-rmax, rmax)
    ax.set_xlim(xmin - pad, xmax + 1.15 * pad)

    if annotate:
        lines = [f"{k}: {v}" for k, v in annotate.items()]
        annotation_box(ax, "\n".join(lines), fontsize=8.5)

    ax.legend(loc="upper right", fontsize=8, ncol=2, frameon=True, fancybox=True)
    finalize_figure(fig)
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig
