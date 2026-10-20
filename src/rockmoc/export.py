"""Contour export: CSV (data / CAD point curve), DXF, and OpenFOAM blockMeshDict.

* ``export_csv``        — x, r, θ, M with header (data exchange)
* ``export_cad_points`` — header-less ``x y z`` text for SolidWorks "Curve
  Through XYZ Points" / CATIA GSD point import (z = 0, r → y)
* ``export_dxf``        — ASCII DXF R12 POLYLINE (no dependency), optional
  mirrored half and centreline
* ``export_blockmesh``  — axisymmetric wedge ``system/blockMeshDict``
  (OpenFOAM v2212+/v2306), spline wall edges, collapsed axis patch
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from rockmoc._version import __version__
from rockmoc.geometry.contour import Contour

_UNIT_SCALE = {"m": 1.0, "mm": 1000.0, "in": 1.0 / 0.0254}
_DXF_INSUNITS = {"m": 6, "mm": 4, "in": 1}


def _scale(units: str) -> float:
    if units not in _UNIT_SCALE:
        raise ValueError(f"units must be one of {sorted(_UNIT_SCALE)}")
    return _UNIT_SCALE[units]


def export_csv(contour: Contour, path: str | Path, *, units: str = "m") -> Path:
    """CSV with header ``x,r,theta_rad[,M]`` in the given length units."""
    k = _scale(units)
    cols = [contour.x * k, contour.r * k, contour.theta]
    header = f"x_{units},r_{units},theta_rad"
    if contour.M is not None:
        cols.append(contour.M)
        header += ",M"
    p = Path(path)
    np.savetxt(p, np.column_stack(cols), delimiter=",", header=header, comments="", fmt="%.9g")
    return p


def export_cad_points(contour: Contour, path: str | Path, *, units: str = "mm") -> Path:
    """Header-less ``x y z`` point list (z = 0) for CAD curve-through-points import."""
    k = _scale(units)
    pts = np.column_stack([contour.x * k, contour.r * k, np.zeros(len(contour.x))])
    p = Path(path)
    np.savetxt(p, pts, delimiter=" ", fmt="%.9g")
    return p


def _dxf_polyline(xs: np.ndarray, ys: np.ndarray, layer: str) -> list[str]:
    out = ["0", "POLYLINE", "8", layer, "66", "1", "70", "0"]
    for x, y in zip(xs, ys):
        out += ["0", "VERTEX", "8", layer, "10", f"{x:.9g}", "20", f"{y:.9g}", "30", "0.0"]
    out += ["0", "SEQEND", "8", layer]
    return out


def export_dxf(
    contour: Contour,
    path: str | Path,
    *,
    units: str = "mm",
    mirror: bool = False,
    centerline: bool = True,
    n_points: int | None = None,
) -> Path:
    """ASCII DXF (R12) with the wall as a POLYLINE on layer ``WALL``.

    ``mirror`` adds the r < 0 half; ``centerline`` adds the axis on layer
    ``AXIS``. ``n_points`` resamples by arc length (default: as stored).
    """
    k = _scale(units)
    c = contour.resample(n_points) if n_points else contour
    x, r = c.x * k, c.r * k
    lines = ["999", f"rockmoc {__version__} nozzle contour",
             "0", "SECTION", "2", "HEADER",
             "9", "$ACADVER", "1", "AC1009",
             "9", "$INSUNITS", "70", str(_DXF_INSUNITS[units]),
             "0", "ENDSEC",
             "0", "SECTION", "2", "TABLES",
             "0", "TABLE", "2", "LAYER", "70", "2",
             "0", "LAYER", "2", "WALL", "70", "0", "62", "5", "6", "CONTINUOUS",
             "0", "LAYER", "2", "AXIS", "70", "0", "62", "8", "6", "CONTINUOUS",
             "0", "ENDTAB", "0", "ENDSEC",
             "0", "SECTION", "2", "ENTITIES"]
    lines += _dxf_polyline(x, r, "WALL")
    if mirror:
        lines += _dxf_polyline(x, -r, "WALL")
    if centerline:
        lines += ["0", "LINE", "8", "AXIS", "10", f"{x[0]:.9g}", "20", "0.0", "30", "0.0",
                  "11", f"{x[-1]:.9g}", "21", "0.0", "31", "0.0"]
    lines += ["0", "ENDSEC", "0", "EOF"]
    p = Path(path)
    p.write_text("\n".join(lines) + "\n", encoding="ascii", errors="replace")
    return p


def export_blockmesh(
    contour: Contour,
    path: str | Path,
    *,
    n_axial: int = 200,
    n_radial: int = 60,
    wedge_deg: float = 5.0,
    wall_grading: float = 0.1,
    n_spline: int = 120,
) -> Path:
    """Axisymmetric wedge ``blockMeshDict`` for the divergent section (throat → exit).

    One hex block with the axis collapsed (``axis`` patch, type empty),
    ``front``/``back`` wedge patches, ``inlet`` (throat plane), ``outlet``
    (exit plane) and ``wall`` (spline edges through the contour).
    ``wall_grading`` = last/first radial cell size (< 1 clusters at the wall
    for y⁺ control). Lengths in metres (``scale 1``).
    """
    if not 0.0 < wedge_deg < 10.0:
        raise ValueError("wedge_deg should be in (0, 10) degrees")
    h = math.radians(0.5 * wedge_deg)
    ch, sh = math.cos(h), math.sin(h)
    x0, x1 = float(contour.x[0]), float(contour.x[-1])
    r0, r1 = float(contour.r[0]), float(contour.r[-1])
    c = contour.resample(n_spline + 2)
    inner = list(zip(c.x[1:-1], c.r[1:-1]))

    def v(x: float, y: float, z: float) -> str:
        return f"    ({x:.9g} {y:.9g} {z:.9g})"

    def spline(a: int, b: int, sgn: float) -> str:
        pts = "\n".join(v(x, r * ch, sgn * r * sh) for x, r in inner)
        return f"    spline {a} {b}\n    (\n{pts}\n    )"

    text = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  rockmoc {__version__} — axisymmetric nozzle wedge (divergent section)
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

scale 1;

vertices
(
{v(x0, 0.0, 0.0)}    // 0 axis, throat
{v(x1, 0.0, 0.0)}    // 1 axis, exit
{v(x1, r1 * ch, -r1 * sh)}    // 2 wall, exit, back
{v(x0, r0 * ch, -r0 * sh)}    // 3 wall, throat, back
{v(x0, r0 * ch, r0 * sh)}    // 4 wall, throat, front
{v(x1, r1 * ch, r1 * sh)}    // 5 wall, exit, front
);

blocks
(
    hex (0 1 2 3 0 1 5 4) ({int(n_axial)} {int(n_radial)} 1) simpleGrading (1 {wall_grading:g} 1)
);

edges
(
{spline(3, 2, -1.0)}
{spline(4, 5, 1.0)}
);

boundary
(
    inlet  {{ type patch;  faces ((0 0 4 3)); }}
    outlet {{ type patch;  faces ((1 2 5 1)); }}
    wall   {{ type wall;   faces ((3 4 5 2)); }}
    front  {{ type wedge;  faces ((0 1 5 4)); }}
    back   {{ type wedge;  faces ((0 3 2 1)); }}
    axis   {{ type empty;  faces ((0 1 1 0)); }}
);

mergePatchPairs
(
);

// ************************************************************************* //
"""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    return p
