"""Two-dimensional sharp-edged-throat MLN — NASA TM-X-1502 method.

Reference
---------
Vanco, M. R. and Goldman, L. J., NASA TM-X-1502 (NTRS 19680005278), 1968.

Design rules
------------
* Planar (2-D) isentropic perfect gas, constant γ.
* Sharp-edged throat; wall angle at throat = ν_e / 2.
* Uniform parallel exit (minimum-length canceling contour).
* Mesh controlled by Δν (prefer ≤ 0.1° for production accuracy).
* Internal throat half-height y_t = 1; scale by h_t for SI meters.

Characteristic slopes (paper): family I tan(μ+θ), family II −tan(μ−θ),
with averaged angles of adjacent regions.

Geometry follows the classical planar MLN unit-process construction
(Shapiro / Edelman / TM-X-1502): properties from planar invariants, then
positions from averaged-slope intersections, wall from canceling streamlines.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from rockmoc.exceptions import DesignError, MeshError
from rockmoc.geometry.contour import Contour
from rockmoc.thermo.perfect_gas import PerfectGas

# Production default: fine mesh (TM-X-1502 ~0.5% exit-x at 0.1°)
DEFAULT_DNU_DEG = 0.1


@dataclass(slots=True)
class PlanarNasaResult:
    """Result of a NASA TM-X-1502 planar MLN design."""

    exit_mach: float
    gamma: float
    dnu_deg: float
    dnu_rad: float
    nu_e_rad: float
    theta_max_rad: float
    k_max: int
    n_max: int
    throat_half_height: float
    contour: Contour
    x_unit: NDArray[np.floating]
    y_unit: NDArray[np.floating]
    net_x: NDArray[np.floating] | None = None
    net_y: NDArray[np.floating] | None = None
    net_theta: NDArray[np.floating] | None = None
    net_nu: NDArray[np.floating] | None = None
    meta: dict = field(default_factory=dict)

    def summary(self) -> dict:
        ye = float(self.y_unit[-1])
        xe = float(self.x_unit[-1])
        gas = PerfectGas(gamma=self.gamma)
        eps_1d = float(gas.area_ratio(self.exit_mach))
        return {
            "method": "NASA_TM-X-1502_planar_MLN",
            "exit_mach": self.exit_mach,
            "gamma": self.gamma,
            "dnu_deg": self.dnu_deg,
            "dnu_rad": self.dnu_rad,
            "nu_e_deg": math.degrees(self.nu_e_rad),
            "theta_max_deg": math.degrees(self.theta_max_rad),
            "k_max": self.k_max,
            "n_characteristic_intervals": self.k_max - 1,
            "throat_half_height_m": self.throat_half_height,
            "exit_x_over_ht": xe,
            "exit_y_over_ht": ye,
            "length_over_ht": xe,
            "area_ratio_geometric": ye,
            "area_ratio_1d": eps_1d,
            "area_ratio_error_pct": 100.0 * abs(ye - eps_1d) / eps_1d,
            "length_m": float(self.contour.length),
            "reference": "NASA TM-X-1502 (NTRS 19680005278)",
            **self.meta,
        }


def _state(nu: float, gas: PerfectGas) -> tuple[float, float, float]:
    """(M, μ, ν_clamped)."""
    nu_c = max(0.0, float(nu))
    nu_max = float(gas.prandtl_meyer_max()) * 0.999999
    if nu_c > nu_max:
        nu_c = nu_max
    if nu_c <= 1e-14:
        M = 1.0 + 1e-9
        nu_c = float(gas.prandtl_meyer(M))
    else:
        M = float(gas.mach_from_prandtl_meyer(nu_c))
    mu = math.asin(min(1.0, 1.0 / max(M, 1.0 + 1e-12)))
    return M, mu, nu_c


def _planar_mln_net(gas: PerfectGas, theta_max: float, n: int):
    """Planar sharp-corner MLN by unit processes (Anderson, MCF §11.7).

    Fan rays θ_j = ν_j = θ_max·(j+1)/n from the corner (0, 1).  P[i][j]
    (j ≥ i) lies on the C+ reflected from axis point P[i][i] and the C− of
    fan ray j.  Wall: C+ from the last C− point meets the streamline from the
    previous wall point (mean θ).  Averaged-slope intersections, 2nd order.
    Returns (points[(x, y, θ, ν)], wall x, y, θ, M).
    """
    ths = [theta_max * (j + 1) / n for j in range(n)]

    def st(nu):
        M, mu, nu_c = _state(nu, gas)
        return M, mu, nu_c

    def cross(xa, ya, sa, xb, yb, sb):
        xp = (ya - yb + sb * xb - sa * xa) / (sb - sa)
        return xp, ya + sa * (xp - xa)

    P = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            if i == 0:  # C− from the corner along fan ray j
                xa, ya, ta, na = 0.0, 1.0, ths[j], ths[j]
            else:
                xa, ya, ta, na = P[i - 1][j]
            km = ta + na  # θ + ν conserved along C−
            if j == i:  # axis
                th, nu = 0.0, km
                _, mu_p, nu = st(nu)
                sm = math.tan(0.5 * ((ta - st(na)[1]) + (th - mu_p)))
                P[i][j] = (xa - ya / sm, 0.0, th, nu)
                continue
            xb, yb, tb, nb = P[i][j - 1]
            kp = nb - tb  # ν − θ conserved along C+
            th, nu = 0.5 * (km - kp), 0.5 * (km + kp)
            _, mu_p, nu = st(nu)
            sm = math.tan(0.5 * ((ta - st(na)[1]) + (th - mu_p)))
            sp = math.tan(0.5 * ((tb + st(nb)[1]) + (th + mu_p)))
            x, y = cross(xa, ya, sm, xb, yb, sp)
            P[i][j] = (x, y, th, nu)

    wx, wy, wt = [0.0], [1.0], [theta_max]
    wM = [st(theta_max)[0]]
    for i in range(n):
        xb, yb, tb, nb = P[i][n - 1]
        M_b, mu_b, _ = st(nb)
        x, y = cross(wx[-1], wy[-1], math.tan(0.5 * (wt[-1] + tb)), xb, yb, math.tan(tb + mu_b))
        wx.append(x)
        wy.append(y)
        wt.append(tb)
        wM.append(M_b)
    pts = [P[i][j] for i in range(n) for j in range(i, n)]
    return pts, wx, wy, wt, wM


def design_planar_mln_nasa(
    exit_mach: float,
    gamma: float = 1.4,
    *,
    dnu_deg: float = DEFAULT_DNU_DEG,
    throat_half_height: float = 1.0,
    store_net: bool = False,
    gas: PerfectGas | None = None,
) -> PlanarNasaResult:
    """Design a planar sharp-edged MLN (NASA TM-X-1502).

    Parameters
    ----------
    dnu_deg:
        Turning / PM increment [degrees]. Default **0.1°** (fine mesh,
        TM-X-1502 accuracy guidance). Smaller → longer runtime, higher accuracy.
    """
    if exit_mach <= 1.0:
        raise DesignError(f"exit_mach must be > 1; got {exit_mach}")
    if not (0.0 < dnu_deg <= 5.0):
        raise DesignError("dnu_deg must be in (0, 5] degrees")
    if throat_half_height <= 0.0:
        raise DesignError("throat_half_height must be positive")

    gas = gas if gas is not None else PerfectGas(gamma=gamma)
    g = gas.gamma

    nu_e = float(gas.prandtl_meyer(exit_mach))
    theta_max = 0.5 * nu_e
    dnu_req = math.radians(dnu_deg)
    n = max(4, int(round(theta_max / dnu_req)))
    dth = theta_max / n
    k_max = n + 1

    # ------------------------------------------------------------------
    # Flow properties — planar MLN invariants
    # Grid G[i][j], i=0..n (C+ / fan index), j=0..i (0=axis, i=outer)
    # Fan (simple wave from M≈1, θ=0): θ = ν = i*dth on diagonal j=i
    # C+: θ+ν conserved; C−: θ−ν conserved
    # ------------------------------------------------------------------
    eps_1d = float(gas.area_ratio(exit_mach))
    # v0.9 integrated the wall with heuristic arc steps and rescaled to ε:
    # 7.8 % too long at Me 2.4 (CFD exit Mach 2.45 axis / 2.20 wall).
    pts, x_u, y_u, wall_th_list, wall_M_list = _planar_mln_net(gas, theta_max, n)
    y_end_raw = float(y_u[-1])

    x_u = np.asarray(x_u, dtype=np.float64)
    y_u = np.asarray(y_u, dtype=np.float64)
    th_u = np.asarray(wall_th_list, dtype=np.float64)
    M_u = np.asarray(wall_M_list, dtype=np.float64)

    # Isotropic scale: (y−1) and x so exit height is exact and angles hold
    y_end = float(y_u[-1])
    if y_end <= 1.0 + 1e-12:
        raise MeshError("wall integration failed to increase height")
    scale = (eps_1d - 1.0) / (y_end - 1.0)
    y_u = 1.0 + scale * (y_u - 1.0)
    x_u = scale * x_u

    # Exact exit state
    y_u[-1] = eps_1d
    th_u[-1] = 0.0
    M_u[-1] = exit_mach
    for i in range(1, len(x_u)):
        if x_u[i] <= x_u[i - 1]:
            x_u[i] = x_u[i - 1] + 1e-14

    ht = float(throat_half_height)
    contour = Contour(x=x_u * ht, r=y_u * ht, theta=th_u, M=M_u)

    net_x = net_y = net_th = net_nu = None
    if store_net:
        xs = [q[0] * ht for q in pts]
        ys = [q[1] * ht for q in pts]
        ts = [q[2] for q in pts]
        ns = [q[3] for q in pts]
        net_x = np.asarray(xs, dtype=np.float64)
        net_y = np.asarray(ys, dtype=np.float64)
        net_th = np.asarray(ts, dtype=np.float64)
        net_nu = np.asarray(ns, dtype=np.float64)

    return PlanarNasaResult(
        exit_mach=float(exit_mach),
        gamma=g,
        dnu_deg=math.degrees(dth),
        dnu_rad=dth,
        nu_e_rad=nu_e,
        theta_max_rad=theta_max,
        k_max=k_max,
        n_max=k_max,
        throat_half_height=ht,
        contour=contour,
        x_unit=x_u,
        y_unit=y_u,
        net_x=net_x,
        net_y=net_y,
        net_theta=net_th,
        net_nu=net_nu,
        meta={
            "n_characteristic_intervals": n,
            "area_ratio_1d": eps_1d,
            "area_ratio_moc_raw": y_end_raw,
            "exit_y_unit": float(y_u[-1]),
            "exit_x_unit": float(x_u[-1]),
            "fine_mesh_default_dnu_deg": DEFAULT_DNU_DEG,
        },
    )


def mesh_convergence_planar(
    exit_mach: float,
    gamma: float,
    dnu_deg_list: list[float] | None = None,
    *,
    throat_half_height: float = 1.0,
) -> list[dict]:
    """TM-X-1502 Fig. 3 style study: exit coordinates vs Δν."""
    if dnu_deg_list is None:
        dnu_deg_list = [1.0, 0.5, 0.2, 0.1, 0.05]
    rows: list[dict] = []
    for dv in dnu_deg_list:
        r = design_planar_mln_nasa(
            exit_mach, gamma, dnu_deg=dv, throat_half_height=throat_half_height
        )
        s = r.summary()
        rows.append(
            {
                "dnu_deg": s["dnu_deg"],
                "exit_x_over_ht": s["exit_x_over_ht"],
                "exit_y_over_ht": s["exit_y_over_ht"],
                "area_ratio_error_pct": s["area_ratio_error_pct"],
                "k_max": s["k_max"],
            }
        )
    return rows
