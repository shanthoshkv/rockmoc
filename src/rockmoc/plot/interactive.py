"""Interactive MOC field (Plotly, optional: ``pip install rockmoc[viz]``).

Hover any node to read x, r, θ, ν, M; characteristics and the wall are drawn
as lines. Returns a ``plotly.graph_objects.Figure`` (``fig.show()``,
``fig.write_html(path)``).
"""

from __future__ import annotations

import math
from typing import Any

from rockmoc.geometry.contour import Contour
from rockmoc.moc.types import CharNet, FlowPoint
from rockmoc.plot.theme import DARK, LIGHT


def _require_plotly() -> Any:
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise ImportError("plotly is required: pip install rockmoc[viz]") from exc
    return go


def _line_trace(go: Any, lines: list[list[FlowPoint]], color: str, name: str, max_lines: int) -> Any:
    xs: list[float | None] = []
    rs: list[float | None] = []
    step = max(1, len(lines) // max_lines)
    for ln in lines[::step]:
        xs += [p.x for p in ln] + [None]
        rs += [p.r for p in ln] + [None]
    return go.Scatter(x=xs, y=rs, mode="lines", line=dict(color=color, width=1), name=name,
                      hoverinfo="skip")


def plot_moc_field_interactive(
    contour: Contour,
    net: CharNet | None = None,
    *,
    title: str = "MOC field",
    theme: str = "dark",
    max_lines: int = 80,
) -> Any:
    go = _require_plotly()
    c = DARK if theme == "dark" else LIGHT
    fig = go.Figure()
    if net is not None:
        fig.add_trace(_line_trace(go, net.interior, c["c_plus"], "C+", max_lines))
        fig.add_trace(_line_trace(go, net.c_minus, c["c_minus"], "C−", max_lines))
        pts = list({id(p): p for p in [*net.points, *net.wall, *net.axis]}.values())
        fig.add_trace(go.Scatter(
            x=[p.x for p in pts], y=[p.r for p in pts], mode="markers", name="nodes",
            marker=dict(size=4, color=[p.M for p in pts], colorscale="Turbo",
                        colorbar=dict(title="M")),
            customdata=[(math.degrees(p.theta), math.degrees(p.nu), p.M, p.kind.name) for p in pts],
            hovertemplate=("x=%{x:.5g} m<br>r=%{y:.5g} m<br>θ=%{customdata[0]:.3f}°"
                           "<br>ν=%{customdata[1]:.3f}°<br>M=%{customdata[2]:.4f}"
                           "<br>%{customdata[3]}<extra></extra>"),
        ))
    fig.add_trace(go.Scatter(x=list(contour.x), y=list(contour.r), mode="lines", name="wall",
                             line=dict(color=c["wall"], width=3)))
    fig.update_layout(
        title=title, template="plotly_dark" if theme == "dark" else "plotly_white",
        xaxis_title="x [m]", yaxis_title="r [m]", yaxis=dict(scaleanchor="x", scaleratio=1),
    )
    return fig
