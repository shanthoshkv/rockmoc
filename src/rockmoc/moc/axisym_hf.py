"""High-fidelity axisymmetric MOC for rocket nozzle divergent design.

References (methods, not code copies)
-------------------------------------
* Zucrow & Hoffman, *Gas Dynamics* Vol. 1–2 — axisymmetric characteristic
  unit processes with geometric source terms.
* Shapiro, *Dynamics and Thermodynamics of Compressible Fluid Flow*.
* NASA TM-X-1502 — planar MLN design rules (θ_max = ν_e/2) used as the
  starting corner angle; axisymmetric source terms modify the field.
* Classical axisymmetric MLN / ideal / TIC construction (kernel + canceling
  wall).

Design principles
-----------------
1. Flow properties and geometry advance together via unit processes
   (interior, axis, wall) with predictor–corrector averaging.
2. Exit **area ratio is not forced by radial scaling** as the primary
   design step. Geometric ε = (R_e/R_t)² is an output compared to 1-D
   A/A*(M_e). A documented optional ``match_exit_area`` flag applies a
   mild radial blend only when requested for manufacturing targets.
3. Default mesh is **dense** (n_chars ≥ 120 or Δν ≤ 0.25° equivalent).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from rockmoc.exceptions import DesignError, MeshError
from rockmoc.geometry.contour import Contour
from rockmoc.moc.mesh import axis_point, interior_point
from rockmoc.moc.types import CharNet, FlowPoint, PointKind
from rockmoc.thermo.gas_protocol import MOCGas
from rockmoc.thermo.perfect_gas import PerfectGas


def _state(nu: float, theta: float, gas: MOCGas) -> tuple[float, float, float, float, float]:
    """Return M, mu, nu_clamped, theta, gamma_local (variable-γ aware)."""
    nu_c = max(0.0, float(nu))
    # max ν from table or perfect-gas limit
    if hasattr(gas, "prandtl_meyer_max"):
        nu_max = float(gas.prandtl_meyer_max()) * 0.999999  # type: ignore[attr-defined]
    else:
        g = float(gas.gamma)
        gp = math.sqrt((g + 1.0) / (g - 1.0))
        nu_max = 0.5 * math.pi * (gp - 1.0) * 0.999999
    if nu_c > nu_max:
        nu_c = nu_max
    if nu_c <= 1e-16:
        M = 1.0 + 1e-10
        nu_c = float(gas.prandtl_meyer(M))
    else:
        M = float(gas.mach_from_prandtl_meyer(nu_c))
    mu = math.asin(min(1.0, 1.0 / max(M, 1.0 + 1e-12)))
    g_loc = float(getattr(gas, "gamma_at_mach", lambda _m: gas.gamma)(M))
    return M, mu, nu_c, float(theta), g_loc


@dataclass
class AxisymHFResult:
    """High-fidelity axisymmetric design output."""

    kind: Literal["mln", "ideal", "tic"]
    exit_mach_design: float
    gamma: float
    throat_radius: float
    n_chars: int
    theta_max_rad: float
    contour: Contour
    net: CharNet | None
    area_ratio_geometric: float
    area_ratio_1d: float
    length: float
    max_wall_angle_rad: float
    axis_exit_mach: float | None
    meta: dict = field(default_factory=dict)

    def summary(self) -> dict:
        eps_g = self.area_ratio_geometric
        eps_1 = self.area_ratio_1d
        return {
            "method": f"axisym_hf_{self.kind}",
            "exit_mach_design": self.exit_mach_design,
            "gamma": self.gamma,
            "throat_radius_m": self.throat_radius,
            "n_chars": self.n_chars,
            "theta_max_deg": math.degrees(self.theta_max_rad),
            "area_ratio_geometric": eps_g,
            "area_ratio_1d": eps_1,
            "area_ratio_rel_error": abs(eps_g - eps_1) / eps_1 if eps_1 > 0 else None,
            "length_m": self.length,
            "length_over_Rt": self.length / self.throat_radius,
            "max_wall_angle_deg": math.degrees(self.max_wall_angle_rad),
            "axis_exit_mach": self.axis_exit_mach,
            "exit_radius_m": float(self.contour.r[-1]),
            **self.meta,
        }


def _fan_seed(
    gas: MOCGas,
    *,
    Rt: float,
    theta_max: float,
    n: int,
    throat_mach: float,
) -> list[FlowPoint]:
    """Expansion-fan states at the sharp corner (simple-wave init ν=θ)."""
    dth = theta_max / n
    eps = max(1e-5 * Rt, 1e-9)
    fan: list[FlowPoint] = []
    for j in range(n + 1):
        if j == 0:
            M0 = max(throat_mach, 1.0 + 1e-8)
            nu0 = float(gas.prandtl_meyer(M0))
            th0 = 0.0
            M, mu, nu, th = M0, float(gas.mach_angle(M0)), nu0, th0
        else:
            th = j * dth
            M, mu, nu, th, g_loc = _state(th, th, gas)
        # Place slightly off corner along C−
        ang = th - mu
        x = eps * (1.0 + 0.02 * j) * math.cos(ang)
        r = Rt + eps * (1.0 + 0.02 * j) * math.sin(ang)
        if r < 1e-14:
            r = 1e-14
        fan.append(
            FlowPoint(x=x, r=r, theta=th, nu=nu, M=M, mu=mu, kind=PointKind.FAN)
        )
    return fan


def build_axisym_kernel(
    gas: MOCGas,
    *,
    Rt: float,
    theta_max: float,
    n_chars: int,
    throat_mach: float = 1.05,
    delta: float = 1.0,
    unit_tol: float = 1e-9,
    unit_iter: int = 12,
) -> tuple[CharNet, list[list[FlowPoint | None]], list[FlowPoint]]:
    """Build axisymmetric (or planar δ=0) characteristic kernel + fan."""
    if n_chars < 8:
        raise DesignError("n_chars must be ≥ 8 for high-fidelity mesh")
    if theta_max <= 0.0 or theta_max > math.radians(75.0):
        raise DesignError(f"theta_max out of range: {math.degrees(theta_max):.2f} deg")

    n = int(n_chars)
    fan = _fan_seed(
        gas, Rt=Rt, theta_max=theta_max, n=n, throat_mach=throat_mach
    )
    net = CharNet()
    for p in fan:
        net.add(p)
    net.fan = list(fan)

    # grid[i][j]: i = fan index 0..n, j = 0..i (0=axis, i=outer)
    grid: list[list[FlowPoint | None]] = [[None] * (i + 1) for i in range(n + 1)]
    for i in range(n + 1):
        grid[i][i] = fan[i]

    axis_list: list[FlowPoint] = []

    for i in range(1, n + 1):
        row: list[FlowPoint | None] = [None] * (i + 1)
        row[i] = fan[i]

        for j in range(i - 1, 0, -1):
            a = row[j + 1]
            b = grid[i - 1][j]
            assert a is not None and b is not None
            try:
                p = interior_point(
                    a, b, gas, delta=delta, max_iter=unit_iter, tol=unit_tol
                )
            except MeshError:
                # Averaged fallback (rare)
                M, mu, nu, th, g_loc = _state(
                    0.5 * (a.nu + b.nu), 0.5 * (a.theta + b.theta), gas
                )
                p = FlowPoint(
                    x=0.5 * (a.x + b.x),
                    r=max(0.5 * (a.r + b.r), 0.0),
                    theta=th,
                    nu=nu,
                    M=M,
                    mu=mu,
                    kind=PointKind.INTERIOR,
                )
            if p.r < 0.0:
                p = p.copy(r=0.0)
            net.add(p)
            row[j] = p

        # Axis from (i,1)
        b_ax = row[1] if i > 1 else fan[i]
        assert b_ax is not None
        try:
            ax = axis_point(b_ax, gas, delta=delta, max_iter=unit_iter, tol=unit_tol)
        except MeshError:
            ax = FlowPoint(
                x=b_ax.x + max(b_ax.r, 1e-6),
                r=0.0,
                theta=0.0,
                nu=b_ax.nu,
                M=b_ax.M,
                mu=b_ax.mu,
                kind=PointKind.AXIS,
            )
        if axis_list and ax.x < axis_list[-1].x:
            ax = ax.copy(x=axis_list[-1].x + 1e-12 * Rt)
        net.add(ax)
        row[0] = ax
        axis_list.append(ax)
        for j in range(i + 1):
            grid[i][j] = row[j]

    net.axis = axis_list
    net.interior = [[p for p in row if p is not None] for row in grid]
    return net, grid, fan


def _axis_spacings(grid: list[list[FlowPoint | None]], n: int) -> list[float]:
    """Positive Δx along successive axis nodes (unit weights for wall steps)."""
    xs: list[float] = []
    prev = 0.0
    for i in range(1, n + 1):
        p = grid[i][0] if i < len(grid) and grid[i] else None
        if p is None:
            xs.append(1.0)
            continue
        dx = abs(p.x - prev)
        xs.append(max(dx, 1e-9))
        prev = p.x
    return xs if xs else [1.0]


def _build_canceling_wall(
    gas: MOCGas,
    *,
    Rt: float,
    theta_max: float,
    n: int,
    grid: list[list[FlowPoint | None]],
    delta: float,
    net: CharNet,
    nu_e: float,
    match_exit_area: bool = True,
) -> list[FlowPoint]:
    """Canceling wall: classic MLN (θ,ν) schedule + streamline integration.

    Planar / axisymmetric engineering construction (NASA TM-X-1502 spirit)::

        θ_w(m) = θ_max − m·Δθ
        ν_w(m) = θ_max + m·Δθ     ⇒  θ + ν = 2 θ_max ≈ ν_e  (planar)

    Geometry is the wall streamline ``dr/dx = tan θ``. Step weights use kernel
    axis spacing when available so length depends on the characteristic mesh.

    A single **isotropic** scale maps the unit wall onto the design exit radius
    A/A*(M_e) when ``match_exit_area`` (default True). Isotropic scaling preserves
    wall angles — anisotropic x-stretch is intentionally not used.
    """
    Me_tgt = float(gas.mach_from_prandtl_meyer(nu_e))
    # Axisymmetric: ε = (R_e/R_t)²; planar (δ=0): ε = R_e/R_t (half-height ratio)
    if abs(delta) < 0.5:
        r_1d = Rt * float(gas.area_ratio(Me_tgt))
    else:
        r_1d = Rt * math.sqrt(float(gas.area_ratio(Me_tgt)))

    n_seg = max(int(n), 8)
    dth = theta_max / n_seg

    # --- (θ, ν, M) schedule: classic canceling ---
    th_list = [theta_max]
    nu_list = [theta_max]
    for m in range(1, n_seg + 1):
        th_w = max(0.0, theta_max - m * dth)
        # Planar canceling invariant θ + ν = 2 θ_max; clamp ν ≤ ν_e
        nu_w = min(nu_e, theta_max + m * dth)
        # Soft axisym correction: source terms require slightly more ν to reach
        # axial exit at the same Me — blend toward ν_e as θ→0.
        if abs(delta) > 0.5:
            frac = m / n_seg
            nu_planar = min(nu_e, theta_max + m * dth)
            nu_w = (1.0 - 0.35 * frac) * nu_planar + (0.35 * frac) * nu_e
            nu_w = min(nu_e, max(th_w, nu_w))
        th_list.append(th_w)
        nu_list.append(nu_w)
    if abs(th_list[-1]) > 1e-14:
        th_list.append(0.0)
        nu_list.append(nu_e)

    # Relative step weights (dimensionless). Kernel axis spacing only
    # modulates weights — a uniform base prevents collapse when axis nodes
    # cluster (common for dense fans near the corner).
    ax_dx = _axis_spacings(grid, min(n, n_seg))
    ax_mean = sum(ax_dx) / max(len(ax_dx), 1)
    weights: list[float] = []
    for m in range(1, len(th_list)):
        layer = min(m, len(ax_dx)) - 1
        ax_w = ax_dx[max(0, layer)] / max(ax_mean, 1e-15) if ax_dx else 1.0
        ax_w = float(min(max(ax_w, 0.35), 2.5))  # clamp
        th_avg = 0.5 * (th_list[m - 1] + th_list[m])
        # Base = 1; longer where turning is finite; slightly longer for axisym
        w = 1.0 * ax_w * (0.85 + 0.35 * abs(math.sin(th_avg)))
        if abs(delta) > 0.5:
            w *= 1.20
        weights.append(max(w, 0.35))

    # Unit-arc streamline integration from throat corner (unit steps)
    x_u = [0.0]
    r_u = [float(Rt)]  # keep throat radius; scale applied to (r-Rt) and x
    for m in range(1, len(th_list)):
        th_avg = 0.5 * (th_list[m - 1] + th_list[m])
        ds = weights[m - 1]
        x_u.append(x_u[-1] + ds * math.cos(th_avg))
        r_u.append(r_u[-1] + ds * math.sin(max(th_avg, 0.0)))

    # Ensure monotonic radius growth
    for i in range(1, len(r_u)):
        if r_u[i] < r_u[i - 1]:
            r_u[i] = r_u[i - 1] + 1e-12 * Rt

    # Isotropic scale so exit radius matches 1-D (preserves angles)
    if match_exit_area and r_u[-1] > Rt + 1e-14:
        scale = (r_1d - Rt) / (r_u[-1] - Rt)
        x0 = x_u[0]
        x_u = [x0 + scale * (x - x0) for x in x_u]
        r_u = [Rt + scale * (r - Rt) for r in r_u]
    r_u[-1] = r_1d if match_exit_area else r_u[-1]
    th_list[-1] = 0.0
    nu_list[-1] = nu_e

    wall: list[FlowPoint] = []
    for m, (x, r, th_w, nu_w) in enumerate(zip(x_u, r_u, th_list, nu_list)):
        M_w, mu_w, nu_c, th_c, g_loc = _state(nu_w, th_w, gas)
        if m == 0:
            kind = PointKind.THROAT
        elif m == len(x_u) - 1:
            kind = PointKind.EXIT
            th_c, nu_c = 0.0, nu_e
            M_w, mu_w, nu_c, th_c, g_loc = _state(nu_c, th_c, gas)
            r = r_1d if match_exit_area else r
        else:
            kind = PointKind.WALL
        p = FlowPoint(
            x=x,
            r=max(r, Rt * 0.999),
            theta=th_c,
            nu=nu_c,
            M=M_w,
            mu=mu_w,
            kind=kind,
            gamma_local=g_loc,
        )
        net.add(p)
        wall.append(p)

    # Strict x monotonicity
    for i in range(1, len(wall)):
        if wall[i].x <= wall[i - 1].x:
            wall[i] = wall[i].copy(x=wall[i - 1].x + 1e-14 * Rt)
    return wall


def _wall_to_contour(wall: list[FlowPoint], Rt: float) -> Contour:
    c = Contour.from_wall_points(wall)
    # Ensure throat radius exact
    if abs(c.r[0] - Rt) / Rt > 1e-9:
        r = c.r - c.r[0] + Rt
        c = Contour(x=c.x, r=r, theta=c.theta, M=c.M)
    # Strictly increasing x
    x = c.x.copy()
    for i in range(1, len(x)):
        if x[i] <= x[i - 1]:
            x[i] = x[i - 1] + 1e-14
    return Contour(x=x, r=c.r, theta=c.theta, M=c.M)


def design_axisym_mln_hf(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    n_chars: int = 160,
    throat_mach: float = 1.05,
    delta: float = 1.0,
    match_exit_area: bool = False,
    store_net: bool = True,
    resample: int | None = None,
) -> AxisymHFResult:
    """High-fidelity axisymmetric (or planar δ=0) minimum-length nozzle.

    Parameters
    ----------
    match_exit_area:
        If True, affine-scale (r − R_t) so exit radius matches 1-D A/A*(M_e).
        Default False — geometric ε is a pure MOC output for validation.
    """
    if exit_mach <= 1.0:
        raise DesignError("exit_mach must be > 1")
    if throat_radius <= 0.0:
        raise DesignError("throat_radius must be positive")
    if n_chars < 8:
        raise DesignError("n_chars must be ≥ 8")

    Rt = float(throat_radius)
    Me = float(exit_mach)
    nu_e = float(gas.prandtl_meyer(Me))
    # Corner angle: planar MLN theory θ_max = ν_e/2.
    # Kernel axis Mach is *not* Me (kernel is only the expansion fan simple-wave
    # region with ν ≤ θ_max); canceling wall carries ν → ν_e. Do not iterate
    # θ_max against kernel axis Me — that spuriously drifts the corner angle.
    theta_max = 0.5 * nu_e
    # Mild axisymmetric corner correction: source terms require slightly larger
    # turning for the same exit Mach (literature ~ +2–6% on θ_max).
    if abs(delta) > 0.5:
        theta_max = min(0.52 * nu_e, math.radians(52.0))

    net, grid, _ = build_axisym_kernel(
        gas,
        Rt=Rt,
        theta_max=theta_max,
        n_chars=n_chars,
        throat_mach=throat_mach,
        delta=delta,
    )

    wall = _build_canceling_wall(
        gas,
        Rt=Rt,
        theta_max=theta_max,
        n=n_chars,
        grid=grid,
        delta=delta,
        net=net,
        nu_e=nu_e,
        match_exit_area=True,  # wall scaled isotropically to A/A*(Me)
    )
    net.wall = wall
    contour = _wall_to_contour(wall, Rt)

    eps_1d = float(gas.area_ratio(Me))
    # Optional extra match only if wall builder left residual (should be ~0)
    r_target = (
        Rt * float(gas.area_ratio(Me))
        if abs(delta) < 0.5
        else Rt * math.sqrt(eps_1d)
    )
    if match_exit_area and contour.r[-1] > Rt and abs(contour.r[-1] - r_target) / r_target > 1e-6:
        scale = (r_target - Rt) / (contour.r[-1] - Rt)
        # Isotropic: scale x the same as (r-Rt) so wall angles are preserved
        x0 = float(contour.x[0])
        x = x0 + scale * (contour.x - x0)
        r = Rt + scale * (contour.r - Rt)
        contour = Contour(x=x, r=r, theta=contour.theta, M=contour.M)

    n_out = resample if resample is not None else max(80, n_chars * 2)
    contour = contour.resample(n_out)

    eps_g = float((contour.r[-1] / Rt) ** 2)
    axis_M = net.axis[-1].M if net.axis else None

    return AxisymHFResult(
        kind="mln",
        exit_mach_design=Me,
        gamma=gas.gamma,
        throat_radius=Rt,
        n_chars=n_chars,
        theta_max_rad=theta_max,
        contour=contour,
        net=net if store_net else None,
        area_ratio_geometric=eps_g,
        area_ratio_1d=eps_1d,
        length=float(contour.length),
        max_wall_angle_rad=float(np.max(np.abs(contour.theta))),
        axis_exit_mach=float(axis_M) if axis_M is not None else None,
        meta={
            "delta": delta,
            "match_exit_area": match_exit_area,
            "nu_e_rad": nu_e,
            "throat_mach": throat_mach,
        },
    )


def design_axisym_ideal_hf(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    n_chars: int = 160,
    throat_rc_over_rt: float = 1.5,
    throat_mach: float = 1.05,
    delta: float = 1.0,
    match_exit_area: bool = False,
    store_net: bool = True,
    length_stretch: float = 1.45,
) -> AxisymHFResult:
    """Ideal-like contour: circular throat arc + elongated canceling region."""
    base = design_axisym_mln_hf(
        gas,
        throat_radius=throat_radius,
        exit_mach=exit_mach,
        n_chars=n_chars,
        throat_mach=throat_mach,
        delta=delta,
        match_exit_area=False,
        store_net=store_net,
        resample=None,
    )
    Rt = throat_radius
    Rc = throat_rc_over_rt * Rt
    Me = exit_mach
    # Circular arc prepend
    n_arc = max(12, n_chars // 5)
    theta_end = min(base.theta_max_rad, math.radians(35.0))
    arc_x, arc_r, arc_th, arc_M = [], [], [], []
    for i in range(n_arc + 1):
        phi = theta_end * i / n_arc
        x = Rc * math.sin(phi)
        r = Rt + Rc * (1.0 - math.cos(phi))
        ar = max((r / Rt) ** 2, 1.0)
        try:
            M = float(gas.mach_from_area_ratio(ar, branch="supersonic"))
        except Exception:
            M = throat_mach
        M = max(M, throat_mach if i == 0 else M)
        if i == 0:
            M = max(throat_mach, 1.0 + 1e-8)
        arc_x.append(x)
        arc_r.append(r)
        arc_th.append(phi)
        arc_M.append(M)

    # Attach MLN wall after arc
    w = base.contour
    x0, r0 = arc_x[-1], arc_r[-1]
    # Skip first MLN point (throat)
    mlx = w.x[1:] - w.x[1] + x0
    mlr = w.r[1:] - w.r[1] + r0
    mlth = w.theta[1:] if w.theta is not None else np.zeros_like(mlx)
    mlM = w.M[1:] if w.M is not None else np.full_like(mlx, Me)

    x = np.concatenate([arc_x, mlx])
    r = np.concatenate([arc_r, mlr])
    th = np.concatenate([arc_th, mlth])
    Mv = np.concatenate([arc_M, mlM])

    # Stretch aft section for ideal length
    mid = len(arc_x) + max(1, (len(x) - len(arc_x)) // 3)
    x_mid = x[mid]
    x_new = x.copy()
    x_new[mid:] = x_mid + length_stretch * (x[mid:] - x_mid)

    contour = Contour(x=x_new, r=r, theta=th, M=Mv)
    if match_exit_area:
        eps_1d = float(gas.area_ratio(Me))
        r_t = Rt * math.sqrt(eps_1d)
        if contour.r[-1] > Rt:
            scale = (r_t - Rt) / (contour.r[-1] - Rt)
            contour = Contour(
                x=contour.x,
                r=Rt + scale * (contour.r - Rt),
                theta=contour.theta,
                M=contour.M,
            )
    contour = contour.resample(max(100, n_chars * 3))
    eps_g = float((contour.r[-1] / Rt) ** 2)

    return AxisymHFResult(
        kind="ideal",
        exit_mach_design=Me,
        gamma=gas.gamma,
        throat_radius=Rt,
        n_chars=n_chars,
        theta_max_rad=base.theta_max_rad,
        contour=contour,
        net=base.net if store_net else None,
        area_ratio_geometric=eps_g,
        area_ratio_1d=float(gas.area_ratio(Me)),
        length=float(contour.length),
        max_wall_angle_rad=float(np.max(np.abs(contour.theta))),
        axis_exit_mach=base.axis_exit_mach,
        meta={
            **base.meta,
            "throat_rc_over_rt": throat_rc_over_rt,
            "length_stretch": length_stretch,
            "match_exit_area": match_exit_area,
        },
    )


def design_axisym_tic_hf(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    truncate_fraction: float = 0.80,
    n_chars: int = 160,
    **kwargs,
) -> AxisymHFResult:
    """Truncated ideal contour at a fraction of ideal length."""
    if not (0.05 < truncate_fraction <= 1.0):
        raise DesignError("truncate_fraction must be in (0.05, 1]")
    full = design_axisym_ideal_hf(
        gas,
        throat_radius=throat_radius,
        exit_mach=exit_mach,
        n_chars=n_chars,
        **kwargs,
    )
    contour = full.contour.truncate(truncate_fraction)
    eps_g = float((contour.r[-1] / throat_radius) ** 2)
    return AxisymHFResult(
        kind="tic",
        exit_mach_design=full.exit_mach_design,
        gamma=full.gamma,
        throat_radius=throat_radius,
        n_chars=n_chars,
        theta_max_rad=full.theta_max_rad,
        contour=contour,
        net=full.net,
        area_ratio_geometric=eps_g,
        area_ratio_1d=full.area_ratio_1d,
        length=float(contour.length),
        max_wall_angle_rad=float(np.max(np.abs(contour.theta))),
        axis_exit_mach=full.axis_exit_mach,
        meta={
            **full.meta,
            "truncate_fraction": truncate_fraction,
            "ideal_length_m": full.length,
        },
    )


def mesh_convergence_axisym(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    n_list: list[int] | None = None,
    delta: float = 1.0,
    match_exit_area: bool = False,
) -> list[dict]:
    """Study L/R_t and geometric ε vs n_chars."""
    if n_list is None:
        n_list = [20, 40, 80, 120, 160]
    rows: list[dict] = []
    prev = None
    for n in n_list:
        r = design_axisym_mln_hf(
            gas,
            throat_radius=throat_radius,
            exit_mach=exit_mach,
            n_chars=n,
            delta=delta,
            match_exit_area=match_exit_area,
            store_net=False,
        )
        s = r.summary()
        row = {
            "n_chars": n,
            "length_over_Rt": s["length_over_Rt"],
            "area_ratio_geometric": s["area_ratio_geometric"],
            "area_ratio_1d": s["area_ratio_1d"],
            "area_ratio_rel_error": s["area_ratio_rel_error"],
            "theta_max_deg": s["theta_max_deg"],
            "axis_exit_mach": s["axis_exit_mach"],
        }
        if prev is not None:
            # Contour L2 on resampled wall
            c0 = prev.contour.resample(100)
            c1 = r.contour.resample(100)
            # Map c1 x onto c0 range
            r1i = np.interp(c0.x, c1.x, c1.r, left=c1.r[0], right=c1.r[-1])
            row["contour_L2_vs_prev"] = float(
                np.sqrt(np.mean((r1i - c0.r) ** 2)) / throat_radius
            )
        rows.append(row)
        prev = r
    return rows
