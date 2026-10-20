"""Backward-compatible contour plot entry (dark theme). """

from __future__ import annotations

from typing import Any

from rockmoc.geometry.contour import Contour
from rockmoc.moc.types import CharNet
from rockmoc.plot.moc_field import plot_moc_field
from rockmoc.plot.performance_plots import plot_nozzle_contour
from rockmoc.plot.theme import themed


@themed
def plot_contour(
    contour: Contour,
    *,
    net: CharNet | None = None,
    show: bool = True,
    ax: Any = None,
    mirror: bool = True,
    title: str | None = None,
    save_path: str | None = None,
    **kwargs: Any,
) -> Any:
    """Plot wall contour; if ``net`` is given, draw full MOC field."""
    if net is not None:
        return plot_moc_field(
            contour,
            net,
            title=title or "MOC field & nozzle wall",
            mirror=mirror,
            show=show,
            save_path=save_path,
        )
    return plot_nozzle_contour(
        contour,
        title=title or "Nozzle wall contour",
        mirror=mirror,
        show=show,
        save_path=save_path,
    )
