"""Transition region (Goursat problem) and mass-flow wall (Zucrow & Hoffman)."""

from __future__ import annotations

import math

import numpy as np

from rockmoc.moc._common import _find_xy, _fp, _poly_mass
from rockmoc.moc.types import FlowPoint, PointKind
from rockmoc.thermo.gas_protocol import MOCGas


def _transition_wall(
    gas: MOCGas,
    bases: list[list[FlowPoint]],
    *,
    Rt: float,
    delta: float,
    mdot: float | None = None,
    start: FlowPoint | None = None,
    exit_line: list[FlowPoint] | None = None,
    c_minus_out: list[list[FlowPoint]] | None = None,
) -> tuple[list[FlowPoint], list[list[FlowPoint]]]:
    """Transition region + mass-flow wall (Zucrow & Hoffman).

    ``bases`` are C+ polylines ordered from the axis point E (where M = Me)
    outward; each ends on the last expansion C−. The transition region is the
    Goursat problem between that C− and the uniform exit C+ from E (M = Me,
    θ = 0). The wall is the streamline carrying the sonic-throat mass flow,
    located along each C+ (``mdot`` defaults to the uniform sonic throat).
    ``exit_line`` replaces the uniform exit C+ by a given C+ starting at
    ``bases[0][-1]`` (e.g. Rao's control surface); the exit is then not forced
    axial. Returns (wall from upstream to exit, rows).
    """
    from rockmoc.moc.mesh import interior_point

    planar = abs(delta) < 0.5
    E = bases[0][0]
    if planar and start is not None and exit_line is None:
        # Planar: downstream of the last C− the flow is a simple wave, so each
        # C+ is straight with the state of its base point — exact streamline march.
        wall_p: list[FlowPoint] = []
        rows_p: list[list[FlowPoint]] = []
        prev_w = start
        for row in reversed(bases):
            b = row[-1]
            th_w = 0.0 if b is E else b.theta
            x, y = _find_xy((prev_w.x, prev_w.r), (b.x, b.r), 0.5 * (prev_w.theta + th_w), b.theta + b.mu)
            w = _fp(x, y, th_w, b.nu, gas, PointKind.EXIT if b is E else PointKind.WALL)
            wall_p.append(w)
            rows_p.append(row + [w])
            prev_w = w
        return wall_p, rows_p[::-1]
    n_u = max(len(bases), 8)
    if mdot is None:
        mdot = Rt if planar else math.pi * Rt * Rt
    if exit_line is not None:
        prev = list(exit_line)
        lip = prev[-1]
        n_u = len(prev) - 1
    else:
        ar_e = float(gas.area_ratio(E.M))
        r_e = mdot * ar_e if planar else math.sqrt(mdot * ar_e / math.pi)
        t_e = r_e / math.sin(E.mu)
        prev = [
            FlowPoint(x=E.x + t * math.cos(E.mu), r=t * math.sin(E.mu), theta=0.0, nu=E.nu,
                      M=E.M, mu=E.mu, kind=PointKind.INTERIOR, gamma_local=E.gamma_local)
            for t in np.linspace(0.0, t_e, n_u + 1)
        ]
    ext_rows = [prev]
    for row in bases[1:]:
        cur = [row[-1]]
        for i in range(1, n_u + 1):
            cur.append(interior_point(prev[i], cur[-1], gas, delta=delta))
        ext_rows.append(cur)
        prev = cur

    wall: list[FlowPoint] = []
    rows: list[list[FlowPoint]] = []
    inside: list[int] = []
    for row, ext in zip(bases, ext_rows):
        line = row + ext[1:]
        m = _poly_mass(line, gas, planar=planar)
        j = next((j for j in range(1, len(m)) if m[j] >= mdot), len(m) - 1)
        # Quadratic (Lagrange) inverse interpolation in mass when a third point
        # exists: second-order wall placement; linear otherwise.
        quad = (j + 1 < len(m) and m[j + 1] > m[j] > m[j - 1]
                and 0.5 <= (m[j + 1] - m[j]) / (m[j] - m[j - 1]) <= 2.0)
        idx = [j - 1, j, j + 1] if quad else [j - 1, j]
        wts = []
        for k in idx:
            wk = 1.0
            for q in idx:
                if q != k:
                    wk *= (mdot - m[q]) / (m[k] - m[q])
            wts.append(wk)
        x, r_, th, nu = (sum(wk * getattr(line[k], f) for wk, k in zip(wts, idx))
                         for f in ("x", "r", "theta", "nu"))
        w = _fp(x, r_, th, nu, gas, PointKind.WALL)
        wall.append(w)
        rows.append(line[:j] + [w])
        inside.append(j - len(row))  # ext points 1..this index lie inside the wall
    wall.reverse()
    rows.reverse()
    if exit_line is not None:
        wall[-1] = _fp(lip.x, lip.r, lip.theta, lip.nu, gas, PointKind.EXIT)
    else:
        wall[-1] = _fp(wall[-1].x, r_e, 0.0, E.nu, gas, PointKind.EXIT)
    rows[-1][-1] = wall[-1]
    if c_minus_out is not None:
        # C− columns of the transition grid, wall side first (upstream → downstream)
        for i in range(1, n_u + 1):
            col = [ext[i] for k, ext in enumerate(ext_rows) if i <= inside[k]]
            if len(col) >= 2:
                c_minus_out.append(col[::-1])
    return wall, rows
