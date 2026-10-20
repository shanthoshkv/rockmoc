"""Smooth-throat (circular-arc) characteristic lattice for ideal and Rao contours."""

from __future__ import annotations

import math
from typing import Literal

from rockmoc.exceptions import MeshError, ThermoError
from rockmoc.moc._common import _UNIT_ITER, _UNIT_TOL, _fp, _poly_mass
from rockmoc.moc.equations import source_plus
from rockmoc.moc.transition import _transition_wall
from rockmoc.moc.types import CharNet, FlowPoint, PointKind
from rockmoc.thermo.gas_protocol import MOCGas


def build_ideal_net(gas: MOCGas, *, attach_offset: int = 0, **kwargs) -> CharNet:  # type: ignore[no-untyped-def]
    """Single ideal net; see :func:`build_ideal_nets`."""
    return build_ideal_nets(gas, offsets=(attach_offset,), **kwargs)[0]


def build_ideal_nets(
    gas: MOCGas,
    *,
    Rt: float,
    Rc: float,
    nu_target: float | None = None,
    n_chars: int,
    throat_mach: float = 1.01,
    delta: float = 1.0,
    transition: bool = True,
    initial_line: Literal["sauer", "uniform"] = "sauer",
    offsets: tuple[int, ...] = (0,),
    theta_attach: float | None = None,
    lattice: dict | None = None,
) -> list[CharNet]:
    """Ideal (uniform-exit) contour with a circular-arc throat of radius Rc.

    ``initial_line="sauer"`` starts from Sauer's transonic solution on the
    vertical line where the axis reaches throat_mach; ``"uniform"`` starts
    from uniform θ = 0, M = throat_mach on the throat plane (cruder, but lets
    the arc end before the Sauer line for small Rc / low Me).

    Characteristic lattice P(α, β): C+ lines α, C− lines β, marched along
    diagonals s = α + β from the initial line x = 0 (uniform θ = 0,
    M = throat_mach), whose points are I_j = P(m − j, j).

    * interior P(α, β) from C− through P(α−1, β) and C+ through P(α, β−1)
    * wall    P(α, α+m) on the arc x = Rc sinφ, r = Rt + Rc(1 − cosφ):
      C+ ∩ circle, θ_W = φ, ν_W from C+ compatibility
    * axis    P(β+m, β) from the C− through P(β+m−1, β)

    The lattice upstream of any C− does not depend on where the arc ends, so
    the arc is marched until the C− from a wall node reaches the axis at
    ν ≥ ``nu_target``; that node is the attachment (θ_a = its angle) and E is
    the axis end of its C−. Each entry of ``offsets`` attaches that many wall
    nodes earlier (1 ⇒ E just below the target), so callers can interpolate;
    all nets share one characteristic march. The transition region + mass-flow wall follow (:func:`_transition_wall`).
    """
    from rockmoc.moc.mesh import axis_point, interior_point

    m = max(int(n_chars), 6)
    yc = Rt + Rc
    # Initial data line: Sauer's transonic solution on the vertical line where
    # the axis reaches throat_mach (a vertical line stays spacelike near M = 1),
    # or a uniform near-sonic line on the throat plane.
    from rockmoc.moc.transonic import sauer_line

    M0 = max(float(throat_mach), 1.0 + 1e-6)
    if initial_line == "sauer":
        try:
            g_t = float(gas.gamma_at_mach(1.0))  # throat γ
            x_line, _, line = sauer_line(g_t, Rt, Rc, axis_mach=M0, n=m, delta=delta)
        except ValueError as exc:
            raise MeshError(str(exc)) from exc
    else:
        x_line = 0.0
        line = [(Rt * j / m, M0, 0.0) for j in range(m + 1)]
    P: dict[tuple[int, int], FlowPoint] = {}
    for j, (r, M, th) in enumerate(line):
        kind = PointKind.AXIS if j == 0 else (PointKind.THROAT if j == m else PointKind.INTERIOR)
        P[(m - j, j)] = _fp(x_line, r, th, float(gas.prandtl_meyer(M)), gas, kind)

    def _arc_wall(b: FlowPoint) -> FlowPoint:
        al = b.theta + b.mu
        w = None
        for _ in range(_UNIT_ITER):
            dx, dr = math.cos(al), math.sin(al)
            ox, orr = b.x, b.r - yc
            bq = 2.0 * (ox * dx + orr * dr)
            disc = bq * bq - 4.0 * (ox * ox + orr * orr - Rc * Rc)
            if disc < 0.0:
                raise MeshError("C+ misses the throat arc")
            t = 0.5 * (-bq - math.sqrt(disc))
            x, r = b.x + t * dx, b.r + t * dr
            phi = math.atan2(x, yc - r)
            th_p = 0.5 * (b.theta + phi)
            mu_p = b.mu if w is None else 0.5 * (b.mu + w.mu)
            kp = (b.theta - b.nu) + source_plus(th_p, mu_p, 0.5 * (b.r + r), delta=delta) * (x - b.x)
            w_new = _fp(x, r, phi, max(phi - kp, 0.0), gas, PointKind.WALL)
            if w is not None and abs(w_new.nu - w.nu) < _UNIT_TOL:
                return w_new
            w = w_new
            al = 0.5 * (b.theta + b.mu + w.theta + w.mu)
        return w  # type: ignore[return-value]

    def _point(a: int, b: int) -> FlowPoint:
        if b - a == m:
            return _arc_wall(P[(a, b - 1)])
        if a - b == m:
            return axis_point(P[(a - 1, b)], gas, delta=delta)
        return interior_point(P[(a - 1, b)], P[(a, b - 1)], gas, delta=delta)

    # March diagonals s = α + β with the arc continuing. The attachment is the
    # wall node (β = b, α = b − m) whose C− reaches the axis at ν ≥ ν_target,
    # or the first wall node with θ ≥ theta_attach; marching continues until
    # the axis end of that C− exists.
    if (nu_target is None) == (theta_attach is None):
        raise ValueError("give exactly one of nu_target or theta_attach")
    b_end = None
    for s_ in range(m + 1, m + 1 + 400 * m):
        for a in range(max((s_ - m + 1) // 2, 0), (s_ + m) // 2 + 1):
            b = s_ - a
            if abs(a - b) > m or (a, b) in P:
                continue
            try:
                P[(a, b)] = _point(a, b)
            except (MeshError, ThermoError, KeyError):
                continue  # far ahead on the arc; only needed if E is not reached first
            if b_end is not None:
                continue
            if nu_target is not None and a - b == m and b >= m and P[(a, b)].nu >= nu_target:
                b_end = b
            elif theta_attach is not None and b - a == m and b > m and P[(a, b)].theta >= theta_attach:
                b_end = b
        if b_end is not None and (b_end + m, b_end) in P:
            break
        if not any(k[0] + k[1] == s_ for k in P):
            raise MeshError("characteristic net stopped before reaching the target Mach")
    if b_end is None:
        raise MeshError("target Mach not reached on the axis")
    if lattice is not None:  # expose the raw lattice to callers (e.g. Rao TOC)
        lattice.update(P=P, m=m, b_end=b_end)
    b_found = b_end
    P_all = P
    nets: list[CharNet] = []
    for off in offsets:
        b_end = b_found - off
        if b_end <= m:
            raise MeshError("attachment upstream of the initial data line")

        a_E = b_end + m
        a_star = b_end - m
        Pk = {k: v for k, v in P_all.items() if k[1] <= b_end and k[0] <= a_E}
        wall = [Pk[(aw, aw + m)] for aw in range(0, a_star + 1)]

        net = CharNet()
        rows: list[list[FlowPoint]] = []
        for a in sorted({k[0] for k in Pk}):
            rows.append([Pk[k] for k in sorted((k for k in Pk if k[0] == a), key=lambda k: k[1])])
        net.axis = [Pk[(bb + m, bb)] for bb in range(0, b_end + 1) if (bb + m, bb) in Pk]
        net.wall = list(wall)
        net.interior = rows
        for bb in sorted({k[1] for k in Pk}):
            line = [Pk[k] for k in sorted((k for k in Pk if k[1] == bb), key=lambda k: k[0])]
            if len(line) >= 2:
                net.c_minus.append(line)
        if transition:
            bases: list[list[FlowPoint]] = []
            initial = [Pk[(m - j, j)] for j in range(m + 1)]
            for a in range(a_E, a_star, -1):
                line = [Pk[(a, bb)] for bb in range(0, b_end + 1) if (a, bb) in Pk]
                if a <= m:  # C+ starts on the initial line: include the strip below it
                    line = initial[: m - a] + line
                bases.append(line)
            # Mass across the last C− (E → attachment): the wall streamline then
            # passes exactly through the attachment node.
            c_minus = [Pk[(a, b_end)] for a in range(a_E, a_star - 1, -1)]
            mdot = _poly_mass(c_minus, gas, planar=abs(delta) < 0.5)[-1]
            t_wall, t_rows = _transition_wall(gas, bases, Rt=Rt, delta=delta, mdot=mdot, start=wall[-1],
                                              c_minus_out=net.c_minus)
            net.wall = list(wall) + [w for w in t_wall if w.x > wall[-1].x]
            net.interior = rows + t_rows
        net.points = [p for r_ in net.interior for p in r_]
        nets.append(net)
    return nets


def _attach_angle(net: CharNet, Rt: float, Rc: float) -> float:
    """Angle of the last wall node lying on the throat arc."""
    on_arc = [p for p in net.wall if abs(math.hypot(p.x, Rt + Rc - p.r) - Rc) < 1e-9 * Rc]
    return on_arc[-1].theta if on_arc else 0.0
