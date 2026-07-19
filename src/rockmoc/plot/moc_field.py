"""MOC characteristic-field plots (C+, C−, wall, axis) — dark theme."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from rockmoc.geometry.contour import Contour
from rockmoc.moc.types import CharNet, FlowPoint
from rockmoc.plot.theme import (
    COLORS,
    apply_dark_theme,
    maybe_save,
    maybe_show,
    require_matplotlib,
    style_axes,
)


def _structured_rows(net: CharNet) -> list[list[FlowPoint]]:
    """Prefer net.interior rows; else synthesize from points."""
    if net.interior:
        rows = []
        for row in net.interior:
            pts = [p for p in row if p is not None]
            if len(pts) >= 2:
                pts = sorted(pts, key=lambda p: p.r)
                rows.append(pts)
        if rows:
            return rows
    return []


def _plot_c_families(ax: Any, net: CharNet, *, mirror: bool = True, alpha: float = 0.75) -> None:
    """Draw C− / C+ families from structured kernel mesh when available."""
    rows = _structured_rows(net)

    if rows:
        # C− family: within each row (fixed wave index), axis → outer
        for i, row in enumerate(rows):
            xs = [p.x for p in row]
            rs = [p.r for p in row]
            lab = "C− family" if i == 0 else None
            ax.plot(xs, rs, color=COLORS["c_minus"], lw=1.15, alpha=alpha, label=lab, zorder=3)
            if mirror:
                ax.plot(xs, [-r for r in rs], color=COLORS["c_minus"], lw=0.85, alpha=alpha * 0.4, zorder=2)

        # C+ family: connect points of similar index across rows
        max_j = max(len(r) for r in rows)
        for j in range(max_j):
            xs, rs = [], []
            for row in rows:
                if j < len(row):
                    xs.append(row[j].x)
                    rs.append(row[j].r)
            if len(xs) < 2:
                continue
            lab = "C+ family" if j == 0 else None
            ax.plot(xs, rs, color=COLORS["c_plus"], lw=1.15, alpha=alpha, label=lab, zorder=3)
            if mirror:
                ax.plot(xs, [-r for r in rs], color=COLORS["c_plus"], lw=0.85, alpha=alpha * 0.4, zorder=2)
        return

    # Fallback: Riemann-invariant bins on free points
    if not net.points:
        return
    pts = net.points
    ax.plot(
        [p.x for p in pts],
        [p.r for p in pts],
        ".",
        ms=2.5,
        color=COLORS["accent"],
        alpha=0.35,
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
    jp = np.array([p.theta + p.nu for p in pts])
    jm = np.array([p.theta - p.nu for p in pts])
    xs = np.array([p.x for p in pts])
    rs = np.array([p.r for p in pts])
    for family, inv, color, name in (
        ("+", jp, COLORS["c_plus"], "C+ (approx)"),
        ("-", jm, COLORS["c_minus"], "C− (approx)"),
    ):
        bins = np.linspace(inv.min(), inv.max(), 14)
        labeled = False
        for k in range(len(bins) - 1):
            mask = (inv >= bins[k]) & (inv < bins[k + 1])
            if mask.sum() < 2:
                continue
            order = np.argsort(xs[mask] if family == "+" else rs[mask])
            xx = xs[mask][order]
            rr = rs[mask][order]
            lab = name if not labeled else None
            ax.plot(xx, rr, color=color, lw=0.95, alpha=0.55, label=lab, zorder=3)
            if mirror:
                ax.plot(xx, -rr, color=color, lw=0.7, alpha=0.25, zorder=2)
            labeled = True


def _plot_canceling_chars(
    ax: Any,
    wall: list[FlowPoint],
    *,
    mirror: bool = True,
    alpha: float = 0.55,
    n_draw: int = 28,
) -> None:
    """Approximate canceling C− rays from wall points into the core flow.

    Each wall point emits a left-running characteristic with slope
    tan(θ − μ). Drawn a limited length toward the axis for textbook-style
    canceling-region visualization.
    """
    if len(wall) < 3:
        return
    # subsample wall for clarity
    idx = np.linspace(1, len(wall) - 2, min(n_draw, len(wall) - 2), dtype=int)
    labeled = False
    for i in idx:
        p = wall[int(i)]
        # C− direction: dy/dx = tan(θ − μ)
        beta = p.theta - p.mu
        # avoid near-vertical blow-ups
        if abs(math.cos(beta)) < 1e-3:
            continue
        slope = math.tan(beta)
        # integrate backward a fraction of local radius toward axis
        # go upstream (negative x for typical expanding wall with θ>0, μ>θ)
        L = 0.55 * max(p.r, 1e-6)
        # prefer direction that decreases radius (into the duct)
        # try both signs of length; keep the one that reduces r and stays x-reasonable
        best = None
        for sgn in (-1.0, 1.0):
            dx = sgn * L / math.sqrt(1.0 + slope * slope)
            x1 = p.x + dx
            r1 = p.r + slope * dx
            if r1 < p.r and r1 > 0:
                best = (x1, r1)
                break
        if best is None:
            # fall back: step toward axis in r
            r1 = max(0.05 * p.r, p.r - L * 0.6)
            if abs(slope) > 1e-6:
                x1 = p.x + (r1 - p.r) / slope
            else:
                x1 = p.x - 0.3 * L
            best = (x1, r1)
        x1, r1 = best
        lab = "Canceling C−" if not labeled else None
        ax.plot([p.x, x1], [p.r, r1], color=COLORS["c_minus"], lw=0.85, alpha=alpha, ls="--", label=lab, zorder=3)
        if mirror:
            ax.plot([p.x, x1], [-p.r, -r1], color=COLORS["c_minus"], lw=0.65, alpha=alpha * 0.35, ls="--", zorder=2)
        labeled = True


def _plot_wall_nodes(ax: Any, wall: list[FlowPoint], *, mirror: bool = True) -> None:
    if not wall:
        return
    xs = [p.x for p in wall]
    rs = [p.r for p in wall]
    ax.plot(xs, rs, "o", ms=2.8, color=COLORS["wall"], alpha=0.85, label="Wall nodes", zorder=5)
    if mirror:
        ax.plot(xs, [-r for r in rs], "o", ms=2.0, color=COLORS["wall"], alpha=0.35, zorder=4)


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
    """Dark-themed MOC plot: wall, axis, C+/C− characteristics, fan markers.

    Kernel mesh (C+/C−) is drawn from ``net.interior`` when present.
    Canceling-region C− rays are approximated from wall-state (θ, μ).
    """
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)

    fig, ax = plt.subplots(figsize=figsize)
    style_axes(ax, title=title, xlabel="Axial station  x  [m]", ylabel="Radius  r  [m]")

    # Characteristic net first (under wall)
    if net is not None:
        if net.fan:
            fx = [p.x for p in net.fan]
            fr = [p.r for p in net.fan]
            ax.plot(fx, fr, "o", ms=3.8, color=COLORS["fan"], alpha=0.9, label="Expansion fan", zorder=6)
            if mirror:
                ax.plot(fx, [-r for r in fr], "o", ms=2.8, color=COLORS["fan"], alpha=0.4, zorder=5)
            # fan rays to axis / corner
            corner = net.fan[0]
            for p in net.fan[1:]:
                ax.plot(
                    [corner.x, p.x],
                    [corner.r, p.r],
                    color=COLORS["fan"],
                    lw=0.7,
                    alpha=0.45,
                    zorder=3,
                )
                if mirror:
                    ax.plot(
                        [corner.x, p.x],
                        [-corner.r, -p.r],
                        color=COLORS["fan"],
                        lw=0.55,
                        alpha=0.2,
                        zorder=2,
                    )

        _plot_c_families(ax, net, mirror=mirror)

        if net.axis:
            ax.plot(
                [p.x for p in net.axis],
                [p.r for p in net.axis],
                color=COLORS["axis"],
                lw=1.5,
                ls=":",
                label="Axis nodes",
                zorder=4,
            )

        if net.wall:
            _plot_wall_nodes(ax, net.wall, mirror=mirror)
            if show_canceling:
                _plot_canceling_chars(ax, net.wall, mirror=mirror)

    # Nozzle wall (upper + mirror)
    ax.plot(contour.x, contour.r, color=COLORS["wall"], lw=2.8, label="Nozzle wall", zorder=7)
    if mirror:
        ax.plot(contour.x, -contour.r, color=COLORS["wall"], lw=2.8, zorder=7)
        ax.fill_between(
            contour.x,
            contour.r,
            -contour.r,
            color=COLORS["wall_fill"],
            alpha=0.22,
            zorder=1,
            label="Flow passage",
        )

    # Centerline
    x0, x1 = float(contour.x[0]), float(contour.x[-1])
    ax.axhline(0.0, color=COLORS["muted"], lw=1.0, ls="-", alpha=0.7, label="Centerline", zorder=1)

    # Throat / exit markers
    ax.plot([contour.x[0]], [contour.r[0]], "s", color=COLORS["warn"], ms=8, zorder=8, label="Throat")
    ax.plot([contour.x[-1]], [contour.r[-1]], "D", color=COLORS["good"], ms=8, zorder=8, label="Exit")
    if mirror:
        ax.plot([contour.x[0]], [-contour.r[0]], "s", color=COLORS["warn"], ms=6, zorder=8)
        ax.plot([contour.x[-1]], [-contour.r[-1]], "D", color=COLORS["good"], ms=6, zorder=8)

    ax.set_aspect("equal", adjustable="datalim")
    # padding — include kernel upstream of throat if present
    rmax = float(np.max(np.abs(contour.r))) * 1.15 + 1e-6
    xmin = x0
    if net is not None and net.points:
        xmin = min(xmin, min(p.x for p in net.points))
    pad = 0.04 * (x1 - xmin + 1e-9)
    ax.set_ylim(-rmax, rmax)
    ax.set_xlim(xmin - pad, x1 + 1.2 * pad)

    if annotate:
        lines = [f"{k}: {v}" for k, v in annotate.items()]
        ax.text(
            0.02,
            0.98,
            "\n".join(lines),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            color=COLORS["text"],
            bbox=dict(boxstyle="round,pad=0.4", facecolor=COLORS["bg"], edgecolor=COLORS["grid"], alpha=0.9),
            zorder=10,
        )

    ax.legend(loc="upper right", fontsize=8, ncol=2, frameon=True)
    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig
