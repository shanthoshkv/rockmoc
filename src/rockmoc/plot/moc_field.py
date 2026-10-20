"""MOC characteristic-field plots (C+, C−, wall, axis) — dark theme.

Draws the **full** characteristic net when ``net.interior`` rows span the
nozzle (axis → wall). Fan rays, C− (within-row), C+ (across-row), wall
nodes, and the nozzle contour are layered for textbook-style MOC figures.
"""

from __future__ import annotations

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
    style_colorbar,
    themed,
)


def _draw_lines(ax: Any, lines: list[list[FlowPoint]], *, color: str, label: str,
                mirror: bool, alpha: float, max_lines: int) -> None:
    lines = [ln for ln in lines if len(ln) >= 2]
    if not lines:
        return
    step = max(1, len(lines) // max_lines)
    keep = lines[::step]
    if lines[-1] is not keep[-1]:
        keep.append(lines[-1])
    for i, ln in enumerate(keep):
        xs = [p.x for p in ln]
        rs = [p.r for p in ln]
        ax.plot(xs, rs, color=color, lw=0.9, alpha=alpha, label=label if i == 0 else None, zorder=3)
        if mirror:
            ax.plot(xs, [-r for r in rs], color=color, lw=0.6, alpha=alpha * 0.35, zorder=2)


def _plot_c_families(
    ax: Any,
    net: CharNet,
    *,
    mirror: bool = True,
    alpha: float = 0.85,
    max_lines: int = 60,
) -> None:
    """C+ family = ``net.interior`` rows; C− family = ``net.c_minus`` polylines."""
    _draw_lines(ax, net.interior, color=COLORS["c_plus"], label="C+ characteristics",
                mirror=mirror, alpha=alpha, max_lines=max_lines)
    _draw_lines(ax, net.c_minus, color=COLORS["c_minus"], label="C− characteristics",
                mirror=mirror, alpha=alpha, max_lines=max_lines)
    if not net.interior and not net.c_minus and net.points:
        ax.plot([p.x for p in net.points], [p.r for p in net.points], ".", ms=2.5,
                color=COLORS["accent"], alpha=0.4, label="MOC nodes", zorder=3)


def _plot_mach_field(ax: Any, fig: Any, net: CharNet, contour: Contour, *, mirror: bool) -> None:
    """Filled Mach contours from the net nodes (Delaunay), masked outside the wall."""
    import matplotlib.tri as mtri

    pts = list({id(p): p for p in [*net.points, *net.wall, *net.axis]}.values())
    if len(pts) < 10:
        return
    x = np.array([p.x for p in pts])
    r = np.array([p.r for p in pts])
    M = np.array([p.M for p in pts])
    tri = mtri.Triangulation(x, r)
    xc = x[tri.triangles].mean(axis=1)
    rc = r[tri.triangles].mean(axis=1)
    r_wall = np.interp(xc, contour.x, contour.r, left=np.nan, right=np.nan)
    scale = float(contour.r[-1]) / max(net.wall[-1].r, 1e-12) if net.wall else 1.0
    tri.set_mask(~(rc * scale <= r_wall * 1.001))
    levels = np.linspace(float(M.min()), float(M.max()), 24)
    cf = ax.tricontourf(tri, M, levels=levels, cmap="turbo", alpha=0.55, zorder=1)
    if mirror:
        tri_m = mtri.Triangulation(x, -r, tri.triangles, mask=tri.mask)
        ax.tricontourf(tri_m, M, levels=levels, cmap="turbo", alpha=0.55, zorder=1)
    cbar = fig.colorbar(cf, ax=ax, fraction=0.035, pad=0.02)
    style_colorbar(cbar, "Mach number")


def _plot_regions(ax: Any, net: CharNet, *, mirror: bool) -> None:
    """Shade the uniform-exit region (axis point E → lip → exit plane) when present."""
    if not net.axis or not net.wall or abs(net.wall[-1].theta) > 1e-6:
        return
    E, lip = net.axis[-1], net.wall[-1]
    if lip.x <= E.x:
        return
    for sgn in ((1, -1) if mirror else (1,)):
        ax.fill([E.x, lip.x, lip.x], [0.0, sgn * lip.r, 0.0], color=COLORS["good"], alpha=0.10,
                zorder=1, label="Uniform exit region" if sgn == 1 else None)


def _plot_wall_mach(ax: Any, contour: Contour, *, mirror: bool, lw: float = 3.2) -> None:
    """Wall drawn as a LineCollection coloured by local wall Mach."""
    from matplotlib.collections import LineCollection

    if contour.M is None:
        return
    pts = np.column_stack([contour.x, contour.r])
    segs = np.stack([pts[:-1], pts[1:]], axis=1)
    for sgn in ((1, -1) if mirror else (1,)):
        s_ = segs.copy()
        s_[:, :, 1] *= sgn
        lc = LineCollection(s_, cmap="turbo", linewidths=lw, zorder=7)
        lc.set_array(0.5 * (np.asarray(contour.M)[:-1] + np.asarray(contour.M)[1:]))
        ax.add_collection(lc)


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


@themed
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
    mach_field: bool = True,
    wall_mach: bool = True,
    shade_regions: bool = True,
) -> Any:
    """MOC plot: Mach field, C+/C− characteristics, wall coloured by Mach, regions.

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
        if mach_field:
            _plot_mach_field(ax, fig, net, contour, mirror=mirror)
        if shade_regions:
            _plot_regions(ax, net, mirror=mirror)
        _plot_c_families(ax, net, mirror=mirror)
        if net.fan:  # sharp-corner nets only
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
    if wall_mach and contour.M is not None:
        _plot_wall_mach(ax, contour, mirror=mirror)
        ax.plot([], [], color=COLORS["wall"], lw=2.9, label="Nozzle wall (coloured by M)")
    else:
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

    ax.set_aspect("equal", adjustable="box")
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

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), fontsize=8, ncol=4,
              frameon=True, fancybox=True)
    finalize_figure(fig)
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig
