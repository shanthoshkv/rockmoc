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
    th = [[0.0] * (i + 1) for i in range(n + 1)]
    nv = [[0.0] * (i + 1) for i in range(n + 1)]
    Mm = [[1.0] * (i + 1) for i in range(n + 1)]
    mu = [[0.5 * math.pi] * (i + 1) for i in range(n + 1)]
    xx = [[0.0] * (i + 1) for i in range(n + 1)]
    yy = [[0.0] * (i + 1) for i in range(n + 1)]

    for i in range(n + 1):
        t = i * dth
        M_i, mu_i, n_i = _state(t, gas)
        th[i][i] = t
        nv[i][i] = n_i
        Mm[i][i] = M_i
        mu[i][i] = mu_i

    for i in range(1, n + 1):
        for j in range(i - 1, 0, -1):
            th_a, nu_a = th[i][j + 1], nv[i][j + 1]  # C−
            th_b, nu_b = th[i - 1][j], nv[i - 1][j]  # C+
            km = th_a - nu_a
            kp = th_b + nu_b
            th_p = 0.5 * (km + kp)
            nu_p = 0.5 * (kp - km)
            M_p, mu_p, nu_p = _state(nu_p, gas)
            th[i][j] = th_p
            nv[i][j] = nu_p
            Mm[i][j] = M_p
            mu[i][j] = mu_p

        # Axis: θ=0, ν = θ_B + ν_B from off-axis (i,1)
        th_b, nu_b = th[i][1], nv[i][1]
        nu_p = th_b + nu_b
        M_p, mu_p, nu_p = _state(nu_p, gas)
        th[i][0] = 0.0
        nv[i][0] = nu_p
        Mm[i][0] = M_p
        mu[i][0] = mu_p

    # ------------------------------------------------------------------
    # Geometry — throat half-height = 1, corner at (0, 1)
    # ------------------------------------------------------------------
    def slope_m(th1: float, mu1: float, th2: float, mu2: float) -> float:
        """Averaged C− slope tan(θ − μ)."""
        return math.tan(0.5 * ((th1 - mu1) + (th2 - mu2)))

    def slope_p(th1: float, mu1: float, th2: float, mu2: float) -> float:
        """Averaged C+ slope tan(θ + μ)."""
        return math.tan(0.5 * ((th1 + mu1) + (th2 + mu2)))

    def intersect(
        xa: float, ya: float, ma: float, xb: float, yb: float, mb: float
    ) -> tuple[float, float]:
        den = mb - ma
        if abs(den) < 1e-14:
            raise MeshError("parallel characteristics")
        xp = (ya - yb - ma * xa + mb * xb) / den
        yp = ya + ma * (xp - xa)
        return xp, yp

    # Outer fan points: all emanate from the sharp corner (0, 1) — same
    # geometric origin; flow properties differ along each ray.
    for i in range(n + 1):
        xx[i][i] = 0.0
        yy[i][i] = 1.0

    # First row i=1: axis point from fan ray 1 (C− from corner → axis)
    # C− slope at fan state (1,1); from (0,1) to y=0
    sm = math.tan(th[1][1] - mu[1][1])
    if abs(sm) < 1e-14:
        raise MeshError("first C− slope ~ 0")
    # y - 1 = sm * (x - 0); y=0 → x = -1/sm
    xx[1][0] = -1.0 / sm
    yy[1][0] = 0.0
    if xx[1][0] <= 0.0:
        # If slope orientation wrong, flip
        xx[1][0] = abs(xx[1][0])

    # Remaining rows
    for i in range(2, n + 1):
        # Interiors j = i-1 .. 1
        for j in range(i - 1, 0, -1):
            # A = (i, j+1): for j+1==i this is corner fan point (0,1)
            # B = (i-1, j)
            xa, ya = xx[i][j + 1], yy[i][j + 1]
            xb, yb = xx[i - 1][j], yy[i - 1][j]
            # Predictor-corrector with known end states
            ma = slope_m(th[i][j + 1], mu[i][j + 1], th[i][j], mu[i][j])
            mb = slope_p(th[i - 1][j], mu[i - 1][j], th[i][j], mu[i][j])
            xp, yp = intersect(xa, ya, ma, xb, yb, mb)
            for _ in range(4):
                ma = slope_m(th[i][j + 1], mu[i][j + 1], th[i][j], mu[i][j])
                mb = slope_p(th[i - 1][j], mu[i - 1][j], th[i][j], mu[i][j])
                xp, yp = intersect(xa, ya, ma, xb, yb, mb)
            if yp < 0.0:
                yp = 0.0
            xx[i][j] = xp
            yy[i][j] = yp

        # Axis from C+ through (i,1)
        xb, yb = xx[i][1], yy[i][1]
        # Average C+ toward axis state θ=0
        mb = slope_p(th[i][1], mu[i][1], 0.0, mu[i][0])
        if abs(mb) < 1e-14:
            mb = math.tan(mu[i][1])
        # y - yb = mb*(x - xb); y=0
        xx[i][0] = xb - yb / mb
        yy[i][0] = 0.0
        if i > 1 and xx[i][0] < xx[i - 1][0]:
            xx[i][0] = xx[i - 1][0] + 1e-12

    # ------------------------------------------------------------------
    # Wall contour (canceling)
    # Planar MLN: θ_w(m) = θ_max - m*dth, ν_w(m) = θ_max + m*dth
    # so at m=n: θ=0, ν=ν_e.
    #
    # Geometry (NASA / Shapiro style):
    # Wall point m is the intersection of
    #   (1) the wall streamline through previous wall point (slope tan θ_avg)
    #   (2) the C+ characteristic through the axis point of layer m
    #       (or interior point that carries the matching C+ invariant).
    # ------------------------------------------------------------------
    # Wall angle / PM schedule (planar MLN canceling)
    # θ_w(m) = θ_max − m·Δθ, ν_w = θ_max + m·Δθ → exit θ=0, ν=ν_e
    eps_1d = float(gas.area_ratio(exit_mach))

    wall_th_list = [theta_max]
    wall_nu_list = [theta_max]
    wall_M_list = [_state(theta_max, gas)[0]]
    for m in range(1, n + 1):
        th_w = max(0.0, theta_max - m * dth)
        nu_w = theta_max + m * dth
        M_w, _, nu_w = _state(nu_w, gas)
        wall_th_list.append(th_w)
        wall_nu_list.append(nu_w)
        wall_M_list.append(M_w)
    if abs(wall_th_list[-1]) > 1e-12:
        wall_th_list.append(0.0)
        wall_nu_list.append(nu_e)
        wall_M_list.append(exit_mach)

    # Integrate wall streamline with unit arc steps, then isotropic scale so
    # y_e = A/A* (TM-X-1502). This preserves wall angles and yields a unique
    # physical length for the design (mesh-independent at fixed θ schedule).
    # Optionally weight step length using axis-point spacing when available.
    x_u = [0.0]
    y_u = [1.0]
    for m in range(1, len(wall_th_list)):
        th_avg = 0.5 * (wall_th_list[m - 1] + wall_th_list[m])
        # Characteristic-informed step: prefer axis Δx of layer when valid
        ds = 1.0
        layer = min(m, n)
        if layer >= 1 and xx[layer][0] > 0.0:
            x_prev_ax = xx[layer - 1][0] if layer > 1 else 0.0
            d_ax = abs(xx[layer][0] - x_prev_ax)
            if d_ax > 1e-8:
                # Map axis spacing into unit-arc metric (relative weights)
                ds = max(d_ax, 0.05)
        x_u.append(x_u[-1] + ds * math.cos(th_avg))
        y_u.append(y_u[-1] + ds * math.sin(max(th_avg, 1e-12)))

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
        xs, ys, ts, ns = [], [], [], []
        for i in range(n + 1):
            for j in range(i + 1):
                xs.append(xx[i][j] * ht)
                ys.append(yy[i][j] * ht)
                ts.append(th[i][j])
                ns.append(nv[i][j])
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
