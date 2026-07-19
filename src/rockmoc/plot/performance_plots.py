"""Performance and design-trade graphs — dark theme."""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np

from rockmoc.geometry.contour import Contour
from rockmoc.plot.theme import (
    COLORS,
    apply_dark_theme,
    maybe_save,
    maybe_show,
    require_matplotlib,
    style_axes,
)


def plot_nozzle_contour(
    contour: Contour,
    *,
    title: str = "Nozzle wall contour",
    mirror: bool = True,
    metrics: dict | None = None,
    show: bool = True,
    save_path: str | None = None,
    figsize: tuple[float, float] = (10, 5.5),
) -> Any:
    """Publication-style dark contour plot with optional metric box."""
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    fig, ax = plt.subplots(figsize=figsize)
    style_axes(ax, title=title, xlabel="x [m]", ylabel="r [m]")

    ax.plot(contour.x, contour.r, color=COLORS["wall"], lw=2.8, label="Wall", zorder=4)
    if mirror:
        ax.plot(contour.x, -contour.r, color=COLORS["wall"], lw=2.8, zorder=4)
        ax.fill_between(contour.x, contour.r, -contour.r, color=COLORS["wall_fill"], zorder=1)
    ax.axhline(0, color=COLORS["muted"], lw=0.9, alpha=0.8)
    ax.plot(contour.x[0], contour.r[0], "s", color=COLORS["warn"], ms=8, label="Throat")
    ax.plot(contour.x[-1], contour.r[-1], "D", color=COLORS["good"], ms=8, label="Exit")

    if contour.theta is not None:
        # secondary: wall angle color scatter
        sc = ax.scatter(
            contour.x,
            contour.r,
            c=np.degrees(contour.theta),
            cmap="cool",
            s=12,
            zorder=5,
            alpha=0.85,
        )
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Wall angle [deg]", color=COLORS["text"])
        cbar.ax.yaxis.set_tick_params(color=COLORS["muted"])
        plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color=COLORS["muted"])

    if metrics:
        txt = "\n".join(f"{k}: {v}" for k, v in metrics.items())
        ax.text(
            0.02,
            0.98,
            txt,
            transform=ax.transAxes,
            va="top",
            fontsize=9,
            color=COLORS["text"],
            bbox=dict(boxstyle="round", facecolor=COLORS["bg"], edgecolor=COLORS["grid"], alpha=0.92),
        )

    ax.set_aspect("equal", adjustable="datalim")
    rmax = float(np.max(np.abs(contour.r))) * 1.2 + 1e-9
    ax.set_ylim(-rmax, rmax)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig


def plot_wall_profiles(
    contour: Contour,
    *,
    title: str = "Wall profiles along length",
    show: bool = True,
    save_path: str | None = None,
    figsize: tuple[float, float] = (11, 7),
) -> Any:
    """r(x), θ(x), M(x) subplots."""
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)
    fig.suptitle(title, color=COLORS["text"], fontweight="semibold")

    style_axes(axes[0], ylabel="r [m]")
    axes[0].plot(contour.x, contour.r, color=COLORS["wall"], label="r(x)")
    axes[0].legend(loc="best", fontsize=9)

    style_axes(axes[1], ylabel="θ [deg]")
    if contour.theta is not None:
        axes[1].plot(contour.x, np.degrees(contour.theta), color=COLORS["c_minus"], label="θ(x)")
    axes[1].legend(loc="best", fontsize=9)

    style_axes(axes[2], xlabel="x [m]", ylabel="M [-]")
    if contour.M is not None:
        axes[2].plot(contour.x, contour.M, color=COLORS["c_plus"], label="M(x)")
    axes[2].legend(loc="best", fontsize=9)

    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig


def plot_performance_vs_altitude(
    rows: Sequence[dict],
    *,
    title: str = "Performance vs altitude",
    show: bool = True,
    save_path: str | None = None,
    figsize: tuple[float, float] = (10, 5.5),
) -> Any:
    """Plot Cf and Isp vs altitude from off_design_altitude_sweep rows."""
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    h = np.array([r["altitude_m"] for r in rows], dtype=float) / 1000.0
    cf = np.array([r.get("Cf", r.get("Cf_amb", np.nan)) for r in rows], dtype=float)
    isp = np.array([r.get("Isp_s", r.get("Isp_amb", np.nan)) for r in rows], dtype=float)

    fig, ax1 = plt.subplots(figsize=figsize)
    style_axes(ax1, title=title, xlabel="Altitude [km]", ylabel="Thrust coefficient  $C_F$")
    ax1.plot(h, cf, "o-", color=COLORS["wall"], label="$C_F$", ms=6)
    ax1.tick_params(axis="y", labelcolor=COLORS["wall"])

    ax2 = ax1.twinx()
    ax2.set_facecolor(COLORS["panel"])
    ax2.plot(h, isp, "s--", color=COLORS["c_plus"], label="$I_{sp}$ [s]", ms=6)
    ax2.set_ylabel("$I_{sp}$ [s]", color=COLORS["c_plus"])
    ax2.tick_params(axis="y", colors=COLORS["c_plus"])
    for spine in ax2.spines.values():
        spine.set_color(COLORS["grid"])

    lines1, lab1 = ax1.get_legend_handles_labels()
    lines2, lab2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, lab1 + lab2, loc="best", fontsize=9)
    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig


def plot_of_ratio_sweep(
    rows: Sequence[dict],
    *,
    title: str = "Oxidizer/fuel ratio sweep",
    show: bool = True,
    save_path: str | None = None,
    figsize: tuple[float, float] = (10, 5.5),
) -> Any:
    """Plot Isp / Cf / γ vs OF from of_ratio_sweep rows."""
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    of = np.array([r["OF"] for r in rows], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    fig.suptitle(title, color=COLORS["text"], fontweight="semibold")

    style_axes(axes[0], xlabel="O/F [-]", ylabel="$I_{sp,vac}$ [s]")
    if "Isp_vac_s" in rows[0]:
        axes[0].plot(of, [r["Isp_vac_s"] for r in rows], "o-", color=COLORS["c_plus"])

    style_axes(axes[1], xlabel="O/F [-]", ylabel="$C_{F,vac}$")
    if "Cf_vac" in rows[0]:
        axes[1].plot(of, [r["Cf_vac"] for r in rows], "s-", color=COLORS["wall"])

    style_axes(axes[2], xlabel="O/F [-]", ylabel="γ / ε")
    if "gamma" in rows[0]:
        axes[2].plot(of, [r["gamma"] for r in rows], "^-", color=COLORS["c_minus"], label="γ")
    if "area_ratio" in rows[0]:
        axb = axes[2].twinx()
        axb.plot(of, [r["area_ratio"] for r in rows], "d--", color=COLORS["axis"], label="ε")
        axb.set_ylabel("ε", color=COLORS["axis"])
        axb.tick_params(axis="y", colors=COLORS["axis"])
    axes[2].legend(loc="best", fontsize=8)

    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig


def plot_gamma_sweep(
    rows: Sequence[dict],
    *,
    title: str = "Specific-heat ratio sweep",
    show: bool = True,
    save_path: str | None = None,
) -> Any:
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    g = np.array([r["gamma"] for r in rows], dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    fig.suptitle(title, color=COLORS["text"], fontweight="semibold")
    style_axes(axes[0], xlabel="γ [-]", ylabel="L / $R_t$")
    axes[0].plot(g, [r["L_over_Rt"] for r in rows], "o-", color=COLORS["wall"])
    style_axes(axes[1], xlabel="γ [-]", ylabel="Area ratio ε")
    axes[1].plot(g, [r["area_ratio"] for r in rows], "s-", color=COLORS["c_plus"])
    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig


def plot_area_ratio_sweep(
    rows: Sequence[dict],
    *,
    title: str = "Area-ratio (ε) design sweep",
    show: bool = True,
    save_path: str | None = None,
    figsize: tuple[float, float] = (11, 5),
) -> Any:
    """Isp / Cf / L/Rt vs design area ratio from area_ratio_sweep_rows."""
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    eps = np.array([r.get("area_ratio", r.get("area_ratio_target")) for r in rows], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    fig.suptitle(title, color=COLORS["text"], fontweight="semibold")

    style_axes(axes[0], xlabel="ε [-]", ylabel="$I_{sp,vac}$ [s]")
    if "Isp_vac_s" in rows[0]:
        axes[0].plot(eps, [r["Isp_vac_s"] for r in rows], "o-", color=COLORS["c_plus"], ms=6)

    style_axes(axes[1], xlabel="ε [-]", ylabel="$C_{F,vac}$")
    if "Cf_vac" in rows[0]:
        axes[1].plot(eps, [r["Cf_vac"] for r in rows], "s-", color=COLORS["wall"], ms=6)

    style_axes(axes[2], xlabel="ε [-]", ylabel="L / $R_t$  &  $M_e$")
    if "L_over_Rt" in rows[0]:
        axes[2].plot(eps, [r["L_over_Rt"] for r in rows], "^-", color=COLORS["c_minus"], label="L/Rt", ms=6)
    if "exit_mach" in rows[0]:
        axb = axes[2].twinx()
        axb.plot(eps, [r["exit_mach"] for r in rows], "d--", color=COLORS["axis"], label="$M_e$", ms=5)
        axb.set_ylabel("$M_e$", color=COLORS["axis"])
        axb.tick_params(axis="y", colors=COLORS["axis"])
        for spine in axb.spines.values():
            spine.set_color(COLORS["grid"])
    axes[2].legend(loc="best", fontsize=8)

    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig


def plot_contour_overlay(
    contours: dict[str, Contour],
    *,
    title: str = "Nozzle contour overlay",
    show: bool = True,
    save_path: str | None = None,
    figsize: tuple[float, float] = (10, 5.5),
) -> Any:
    """Overlay multiple wall contours (e.g. MLN vs ideal vs TIC vs Rao)."""
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    fig, ax = plt.subplots(figsize=figsize)
    style_axes(ax, title=title, xlabel="x [m]", ylabel="r [m]")
    palette = [
        COLORS["wall"],
        COLORS["c_plus"],
        COLORS["c_minus"],
        COLORS["axis"],
        COLORS["fan"],
        COLORS["line2"],
    ]
    rmax = 0.0
    for i, (name, c) in enumerate(contours.items()):
        color = palette[i % len(palette)]
        ax.plot(c.x, c.r, color=color, lw=2.2, label=name)
        rmax = max(rmax, float(np.max(np.abs(c.r))))
    ax.axhline(0, color=COLORS["muted"], lw=0.9, alpha=0.8)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_ylim(-0.05 * rmax, rmax * 1.15 + 1e-9)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig


def plot_design_comparison(
    rows: Sequence[dict],
    *,
    title: str = "Design family comparison",
    show: bool = True,
    save_path: str | None = None,
) -> Any:
    """Bar charts comparing MLN / ideal / TIC / Rao metrics."""
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    kinds = [r["kind"] for r in rows]
    x = np.arange(len(kinds))
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5))
    fig.suptitle(title, color=COLORS["text"], fontweight="semibold")
    palette = [COLORS["wall"], COLORS["c_plus"], COLORS["c_minus"], COLORS["axis"], COLORS["fan"]]

    for ax, key, ylab in (
        (axes[0], "L_over_Rt", "L / $R_t$"),
        (axes[1], "area_ratio", "ε"),
        (axes[2], "Cf_vac", "$C_{F,vac}$"),
    ):
        style_axes(ax, ylabel=ylab)
        vals = [r.get(key, 0.0) for r in rows]
        bars = ax.bar(x, vals, color=palette[: len(x)], edgecolor=COLORS["grid"])
        ax.set_xticks(x)
        ax.set_xticklabels(kinds, rotation=15, ha="right")
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height(), f"{v:.3g}", ha="center", va="bottom", fontsize=8, color=COLORS["muted"])

    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig


def plot_mach_pressure_along_wall(
    contour: Contour,
    gas: Any,
    *,
    pc: float = 5e6,
    title: str = "Wall edge Mach & pressure",
    show: bool = True,
    save_path: str | None = None,
) -> Any:
    """M(x) and p/pc(x) along wall."""
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    fig, ax1 = plt.subplots(figsize=(10, 5))
    style_axes(ax1, title=title, xlabel="x [m]", ylabel="Mach number")
    if contour.M is not None:
        ax1.plot(contour.x, contour.M, color=COLORS["c_plus"], label="M(x)")
    ax2 = ax1.twinx()
    if contour.M is not None:
        pe_pc = np.array([float(gas.pe_pc_from_mach(m)) for m in contour.M])
        ax2.plot(contour.x, pe_pc, color=COLORS["c_minus"], ls="--", label="p/p$_c$")
    ax2.set_ylabel("p / p$_c$", color=COLORS["c_minus"])
    ax2.tick_params(axis="y", colors=COLORS["c_minus"])
    lines1, lab1 = ax1.get_legend_handles_labels()
    lines2, lab2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, lab1 + lab2, loc="best", fontsize=9)
    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig


def plot_bl_displacement(
    inviscid: Contour,
    corrected: Contour,
    delta_star: np.ndarray,
    *,
    title: str = "Boundary-layer displacement correction",
    show: bool = True,
    save_path: str | None = None,
) -> Any:
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(title, color=COLORS["text"], fontweight="semibold")
    style_axes(axes[0], xlabel="x [m]", ylabel="r [m]")
    axes[0].plot(inviscid.x, inviscid.r, color=COLORS["muted"], ls="--", label="Inviscid")
    axes[0].plot(corrected.x, corrected.r, color=COLORS["wall"], label="BL-corrected")
    axes[0].legend(fontsize=9)
    axes[0].set_aspect("equal", adjustable="datalim")
    style_axes(axes[1], xlabel="x [m]", ylabel="δ* [m]")
    axes[1].plot(inviscid.x, delta_star, color=COLORS["warn"], label="δ*")
    axes[1].legend(fontsize=9)
    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig


def plot_aerospike(
    X: np.ndarray,
    R: np.ndarray,
    *,
    Re: float | None = None,
    title: str = "Aerospike contour",
    show: bool = True,
    save_path: str | None = None,
) -> Any:
    plt, mpl = require_matplotlib()
    apply_dark_theme(plt, mpl)
    fig, ax = plt.subplots(figsize=(10, 5))
    style_axes(ax, title=title, xlabel="X [m] from lip", ylabel="R [m]")
    ax.plot(X, R, color=COLORS["wall"], lw=2.5, label="Spike")
    ax.plot(X, -R, color=COLORS["wall"], lw=2.5)
    if Re is not None:
        ax.axhline(Re, color=COLORS["fan"], ls=":", label="Shroud lip $R_e$")
        ax.axhline(-Re, color=COLORS["fan"], ls=":")
        ax.plot([0], [Re], "o", color=COLORS["warn"], label="Lip")
    ax.axhline(0, color=COLORS["muted"], lw=0.8)
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(fontsize=9)
    fig.tight_layout()
    maybe_save(fig, save_path)
    maybe_show(fig, show)
    return fig
