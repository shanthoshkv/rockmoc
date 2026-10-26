"""Axisymmetric transonic throat flow by time-marching Euler (Prozan/Saunders approach).

For tight throats (Rc/Rt < 0.4, steep inlets) no local series converges
(Levine & Coats 1970). This solves the steady inviscid flow through the
convergent–throat–early-divergent region directly: cell-centred finite volume
on a wall-fitted grid (r = η r_w(x)), first-order Rusanov fluxes, local time
stepping to steady state, and Richardson extrapolation over two grids
(:func:`solve_throat_extrapolated`). At Rc/Rt = 2 this reproduces the
Kliegel–Levine series (axis M 0.897, Cd 0.9970 vs 0.9962); at 0.625 it gives
axis/wall M 0.82/1.33 and Cd 0.982 against 0.8/1.4 and 0.985 measured
(Cuffel et al. 1969). Boundaries: subsonic inflow with
p0, T0 and grid-aligned flow angle; supersonic extrapolated outflow; slip wall;
symmetry axis (zero-area face). Outputs the throat-plane Mach profile, the
discharge coefficient and a start line for MOC at a chosen x.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np


@dataclass
class EulerThroatResult:
    x: np.ndarray  # cell centres (Nx, Nr)
    r: np.ndarray
    M: np.ndarray
    theta: np.ndarray  # flow angle [rad]
    p_p0: np.ndarray
    Cd: float
    residual: float
    iterations: int

    def plane(self, x0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(r, M, θ) on the column of cells nearest x0."""
        i = int(np.argmin(np.abs(self.x[:, 0] - x0)))
        return self.r[i], self.M[i], self.theta[i]


def solve_throat(
    r_wall: Callable[[np.ndarray], np.ndarray],
    x_in: float,
    x_out: float,
    *,
    gamma: float = 1.4,
    nx: int = 140,
    nr: int = 40,
    cfl: float = 0.5,
    max_iter: int = 30000,
    tol: float = 1e-7,
    x_throat: float = 0.0,
) -> EulerThroatResult:
    """Steady Euler flow in a nozzle r_wall(x), lengths and Rt = r_wall(x_throat) = 1 scale-free.

    Non-dimensional: p0 = 1, ρ0 = 1 (so a0² = γ). Cd = ṁ / (ρ* a* A_t) on the throat plane.
    """
    g = gamma
    xs = np.linspace(x_in, x_out, nx + 1)
    eta = np.linspace(0.0, 1.0, nr + 1)
    rw = r_wall(xs)
    X = np.repeat(xs[:, None], nr + 1, axis=1)
    R = rw[:, None] * eta[None, :]
    # cell geometry
    xc = 0.25 * (X[:-1, :-1] + X[1:, :-1] + X[1:, 1:] + X[:-1, 1:])
    rc = 0.25 * (R[:-1, :-1] + R[1:, :-1] + R[1:, 1:] + R[:-1, 1:])
    area = 0.5 * np.abs((X[1:, 1:] - X[:-1, :-1]) * (R[:-1, 1:] - R[1:, :-1])
                        - (R[1:, 1:] - R[:-1, :-1]) * (X[:-1, 1:] - X[1:, :-1]))
    vol = area * rc
    # i-faces (nx+1, nr): normal (dr, −dx) of the edge (i,j)→(i,j+1); scaled by length and face radius
    ex, er = X[:, 1:] - X[:, :-1], R[:, 1:] - R[:, :-1]
    Li = np.hypot(ex, er)
    nix, nir = er / Li, -ex / Li
    wi = Li * 0.5 * (R[:, 1:] + R[:, :-1])
    # j-faces (nx, nr+1): normal (−dr, dx) of the edge (i,j)→(i+1,j)
    fx, fr = X[1:, :] - X[:-1, :], R[1:, :] - R[:-1, :]
    Lj = np.hypot(fx, fr)
    njx, njr = -fr / Lj, fx / Lj
    wj = Lj * 0.5 * (R[1:, :] + R[:-1, :])
    hmin = np.minimum(np.abs(xs[1:] - xs[:-1])[:, None], (rw[:-1, None] + rw[1:, None]) / (2 * nr))

    # initial guess: 1-D isentropic by area ratio, subsonic upstream of the throat
    from rockmoc.thermo.perfect_gas import PerfectGas

    pg = PerfectGas(g)
    rwc = 0.5 * (rw[:-1] + rw[1:])
    At = float(r_wall(np.array([x_throat]))[0])
    M1 = np.array([float(pg.mach_from_area_ratio(max((rr / At) ** 2, 1.0 + 1e-9),
                                                  branch="subsonic" if xm < x_throat else "supersonic"))
                   for rr, xm in zip(rwc, 0.5 * (xs[:-1] + xs[1:]))])
    T = 1.0 / (1.0 + 0.5 * (g - 1) * M1**2)  # T/T0 ; with ρ0 = p0 = 1, T0 = 1/R → use p = ρ T scaled
    rho = (T ** (1 / (g - 1)))[:, None] * np.ones((1, nr))
    p = (T ** (g / (g - 1)))[:, None] * np.ones((1, nr))
    a = np.sqrt(g * p / rho)
    slope = np.gradient(rw, xs)
    sl_c = 0.5 * (slope[:-1] + slope[1:])[:, None] * (0.5 * (eta[:-1] + eta[1:]))[None, :]
    V = (M1[:, None] * a)
    u, v = V / np.sqrt(1 + sl_c**2), V * sl_c / np.sqrt(1 + sl_c**2)
    U = np.stack([rho, rho * u, rho * v, p / (g - 1) + 0.5 * rho * (u * u + v * v)])

    def prim(Uc):
        rho_ = Uc[0]
        u_, v_ = Uc[1] / rho_, Uc[2] / rho_
        p_ = (g - 1) * (Uc[3] - 0.5 * rho_ * (u_ * u_ + v_ * v_))
        return rho_, u_, v_, p_

    def flux(rho_, u_, v_, p_, nx_, nr_):
        Vn = u_ * nx_ + v_ * nr_
        E = p_ / (g - 1) + 0.5 * rho_ * (u_ * u_ + v_ * v_)
        return np.stack([rho_ * Vn, rho_ * u_ * Vn + p_ * nx_, rho_ * v_ * Vn + p_ * nr_, (E + p_) * Vn]), Vn

    def rusanov(WL, WR, nx_, nr_):
        FL, VnL = flux(*WL, nx_, nr_)
        FR, VnR = flux(*WR, nx_, nr_)
        aL, aR = np.sqrt(g * WL[3] / WL[0]), np.sqrt(g * WR[3] / WR[0])
        s = np.maximum(np.abs(VnL) + aL, np.abs(VnR) + aR)
        EL = WL[3] / (g - 1) + 0.5 * WL[0] * (WL[1] ** 2 + WL[2] ** 2)
        ER = WR[3] / (g - 1) + 0.5 * WR[0] * (WR[1] ** 2 + WR[2] ** 2)
        UL = np.stack([WL[0], WL[0] * WL[1], WL[0] * WL[2], EL])
        UR = np.stack([WR[0], WR[0] * WR[1], WR[0] * WR[2], ER])
        return 0.5 * (FL + FR) - 0.5 * s * (UR - UL)

    def ghosts(W):
        rho_, u_, v_, p_ = W
        # pad i (inflow/outflow) then j (axis/wall); arrays (nx+2, nr+2)
        Wp = [np.pad(q, 1, mode="edge") for q in (rho_, u_, v_, p_)]
        # inflow: p0 = ρ0 = 1, static p extrapolated, direction along the grid line
        pin = np.minimum(p_[0], 0.999)
        Min = np.sqrt(2 / (g - 1) * (pin ** (-(g - 1) / g) - 1))
        Tin = 1 / (1 + 0.5 * (g - 1) * Min**2)
        rin = Tin ** (1 / (g - 1))
        Vin = Min * np.sqrt(g * pin / rin)
        sin_ = slope[0] * 0.5 * (eta[:-1] + eta[1:])
        Wp[0][0, 1:-1], Wp[3][0, 1:-1] = rin, pin
        Wp[1][0, 1:-1], Wp[2][0, 1:-1] = Vin / np.sqrt(1 + sin_**2), Vin * sin_ / np.sqrt(1 + sin_**2)
        # axis: mirror v
        Wp[2][:, 0] = -Wp[2][:, 1]
        # wall: reflect the velocity normal to the wall face
        nxw, nrw = np.pad(njx[:, -1], 1, mode="edge"), np.pad(njr[:, -1], 1, mode="edge")
        uw, vw = Wp[1][:, -2], Wp[2][:, -2]
        vn = uw * nxw + vw * nrw
        Wp[1][:, -1], Wp[2][:, -1] = uw - 2 * vn * nxw, vw - 2 * vn * nrw
        for q in (0, 3):
            Wp[q][:, -1] = Wp[q][:, -2]
        return Wp

    it, res = 0, 1.0
    for it in range(1, max_iter + 1):
        W = prim(U)
        Wp = ghosts(W)
        # first-order Rusanov (second-order MUSCL limit-cycles at the sonic line; accuracy
        # comes from Richardson extrapolation over two grids, see solve_throat_extrapolated)
        Fi = rusanov([q[:-1, 1:-1] for q in Wp], [q[1:, 1:-1] for q in Wp], nix, nir) * wi
        Fj = rusanov([q[1:-1, :-1] for q in Wp], [q[1:-1, 1:] for q in Wp], njx, njr) * wj
        S = np.zeros_like(U)
        S[2] = W[3] * area
        dU = -(Fi[:, 1:] - Fi[:, :-1] + Fj[:, :, 1:] - Fj[:, :, :-1]) + S
        a = np.sqrt(g * W[3] / W[0])
        dt = cfl * hmin / (np.hypot(W[1], W[2]) + a)
        U = U + dt / vol * dU
        if it % 200 == 0:
            res = float(np.max(np.abs(dU[0] * dt / vol / U[0])))
            if res < tol:
                break
    rho_, u_, v_, p_ = prim(U)
    a = np.sqrt(g * p_ / rho_)
    M = np.hypot(u_, v_) / a
    # discharge coefficient on the i-face nearest the throat
    it_ = int(np.argmin(np.abs(xs - x_throat)))
    Wf = [0.5 * (q[it_ - 1] + q[it_]) for q in (rho_, u_, v_, p_)]
    mdot = float(np.sum(Wf[0] * (Wf[1] * nix[it_] + Wf[2] * nir[it_]) * wi[it_]))
    rs, as_ = (2 / (g + 1)) ** (1 / (g - 1)), math.sqrt(g * (2 / (g + 1)))
    Cd = mdot / (rs * as_ * 0.5 * float(rw[it_]) ** 2)
    return EulerThroatResult(xc, rc, M, np.arctan2(v_, u_), p_, Cd, res, it)


def solve_throat_extrapolated(
    r_wall: Callable[[np.ndarray], np.ndarray],
    x_in: float,
    x_out: float,
    *,
    gamma: float = 1.4,
    grids: tuple[tuple[int, int, int], tuple[int, int, int]] = ((160, 44, 9000), (240, 66, 15000)),
    x_throat: float = 0.0,
) -> dict:
    """Throat-plane axis/wall Mach and Cd, Richardson-extrapolated (first-order error ∝ h)."""
    (n1, m1, i1), (n2, m2, i2) = grids
    a = solve_throat(r_wall, x_in, x_out, gamma=gamma, nx=n1, nr=m1, max_iter=i1, x_throat=x_throat)
    b = solve_throat(r_wall, x_in, x_out, gamma=gamma, nx=n2, nr=m2, max_iter=i2, x_throat=x_throat)
    k = 1.0 / (m2 / m1 - 1.0)  # f_exact ≈ f_fine + (f_fine − f_coarse) / (ratio − 1)
    _, Ma, _ = a.plane(x_throat)
    _, Mb, _ = b.plane(x_throat)
    ex = lambda fa, fb: fb + k * (fb - fa)  # noqa: E731
    return {"M_axis": ex(Ma[0], Mb[0]), "M_wall": ex(Ma[-1], Mb[-1]), "Cd": ex(a.Cd, b.Cd), "fine": b}


def cuffel_wall(rc_over_rt: float = 0.625, *, conv_deg: float = 45.0, div_deg: float = 15.0,
                r_inlet: float = 3.0) -> Callable[[np.ndarray], np.ndarray]:
    """Conical convergent–circular-arc throat–conical divergent wall (Rt = 1)."""
    Rc, c, d = rc_over_rt, math.radians(conv_deg), math.radians(div_deg)
    xc, xd = -Rc * math.sin(c), Rc * math.sin(d)
    rcv, rdv = 1 + Rc * (1 - math.cos(c)), 1 + Rc * (1 - math.cos(d))

    def f(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=np.float64)
        arc = 1 + Rc - np.sqrt(np.maximum(Rc**2 - x**2, 0.0))
        out = np.where(x < xc, rcv + (xc - x) * math.tan(c), np.where(x > xd, rdv + (x - xd) * math.tan(d), arc))
        return np.minimum(out, r_inlet)

    return f
