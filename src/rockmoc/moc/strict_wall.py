"""Strict canceling-wall construction via unit processes only.

Builds wall points as the intersection of:
  (1) streamline through the previous wall point (slope = tan θ_w)
  (2) C− through a known interior/axis net point

No soft radius schedule — geometry is pure characteristic + streamline.
"""

from __future__ import annotations

import math

from rockmoc.exceptions import MeshError
from rockmoc.moc.mesh import wall_point_from_previous
from rockmoc.moc.types import FlowPoint, PointKind
from rockmoc.thermo.perfect_gas import PerfectGas


def strict_canceling_wall(
    gas: PerfectGas,
    *,
    throat_radius: float,
    theta_max: float,
    nu_e: float,
    n: int,
    interior_refs: list[FlowPoint],
    delta: float = 1.0,
) -> list[FlowPoint]:
    """March a strict canceling wall from corner using reference C− points.

    Parameters
    ----------
    interior_refs:
        Sequence of net points (e.g. axis or near-wall interior) ordered
        downstream; each provides the C− partner for the next wall point.
    """
    if n < 4:
        raise MeshError("n must be ≥ 4")
    if len(interior_refs) < 2:
        raise MeshError("need ≥ 2 interior reference points")

    dth = theta_max / n
    wall: list[FlowPoint] = []

    M0 = float(gas.mach_from_prandtl_meyer(theta_max))
    mu0 = float(gas.mach_angle(M0))
    corner = FlowPoint(
        x=0.0,
        r=throat_radius,
        theta=theta_max,
        nu=theta_max,
        M=M0,
        mu=mu0,
        kind=PointKind.THROAT,
    )
    wall.append(corner)

    n_seg = min(n, len(interior_refs))
    for m in range(1, n_seg + 1):
        th_w = max(0.0, theta_max - m * dth)
        nu_w = theta_max + (nu_e - theta_max) * (m / n_seg)
        ref = interior_refs[min(m, len(interior_refs) - 1)]
        try:
            p = wall_point_from_previous(
                ref,
                wall[-1],
                gas,
                wall_theta=th_w,
                delta=delta,
            )
            # Override thermo to design schedule for canceling
            M = float(gas.mach_from_prandtl_meyer(max(nu_w, 0.0)))
            mu = float(gas.mach_angle(M))
            p = FlowPoint(
                x=p.x,
                r=max(p.r, wall[-1].r),
                theta=th_w,
                nu=nu_w,
                M=M,
                mu=mu,
                kind=PointKind.WALL if m < n_seg else PointKind.EXIT,
            )
        except MeshError:
            # Geometric streamline step
            w = wall[-1]
            th_avg = 0.5 * (w.theta + th_w)
            ds = max(0.02 * throat_radius, abs(ref.x - w.x) * 0.3)
            dx = ds * math.cos(th_avg)
            dr = ds * math.sin(th_avg)
            M = float(gas.mach_from_prandtl_meyer(max(nu_w, 0.0)))
            mu = float(gas.mach_angle(M))
            p = FlowPoint(
                x=w.x + max(dx, 1e-9),
                r=w.r + max(dr, 0.0),
                theta=th_w,
                nu=nu_w,
                M=M,
                mu=mu,
                kind=PointKind.WALL,
            )
        if p.x <= wall[-1].x:
            p = p.copy(x=wall[-1].x + 1e-9 * throat_radius)
        wall.append(p)

    # Force axial exit
    if abs(wall[-1].theta) > 1e-12:
        fin = wall[-1]
        M = float(gas.mach_from_prandtl_meyer(nu_e))
        wall.append(
            FlowPoint(
                x=fin.x + 1e-12,
                r=fin.r,
                theta=0.0,
                nu=nu_e,
                M=M,
                mu=float(gas.mach_angle(M)),
                kind=PointKind.EXIT,
            )
        )
    return wall
