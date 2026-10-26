"""High-fidelity MOC nozzle design — full characteristic net (kernel + wall).

Algorithm sources (methods, not code copies)
--------------------------------------------
* Anderson, *Modern Compressible Flow* — planar MLN construction
* Zucrow & Hoffman — axisymmetric geometric source terms
* NASA TM-X-1502 — planar sharp-throat rules (θ_max = ν_e/2)
* Planar mesh topology consistent with open references:
  nategphillips/moc, krvi28/MOC_for_Rocket_Nozzle_Design

Design principle
----------------
The **wall is part of the characteristic mesh**. Each wall node inherits the
flow state of the previous interior node (wave cancellation) and is placed at
the intersection of:

1. the C+ characteristic from that interior node, and
2. the wall streamline from the previous wall node.

The full C+/C− net therefore fills the nozzle from throat to exit — the same
mesh used for contour design is drawn in MOC field plots.

Implementation is split across ``_common``, ``kernel_mln``, ``transition``
and ``ideal``; this module holds the public design API and re-exports.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from rockmoc.exceptions import (
    DesignError,
    MeshError,
    MOCFallbackWarning,
    ThermoError,
    warn_fallback,
)
from rockmoc.geometry.contour import Contour
from rockmoc.moc._common import (  # noqa: F401  (re-exported for compatibility)
    _UNIT_ITER,
    _UNIT_TOL,
    _find_xy,
    _fp,
    _isotropic_match_exit,
    _pm_state,
    _poly_mass,
    _records_fallbacks,
    _scale_net,
    _wall_to_contour,
    blend_contours,
)
from rockmoc.moc.ideal import _attach_angle, build_ideal_net, build_ideal_nets  # noqa: F401
from rockmoc.moc.kernel_mln import (  # noqa: F401
    _axisym_transition,
    build_axisym_kernel,
    build_mln_full_net,
)
from rockmoc.moc.transition import _transition_wall  # noqa: F401
from rockmoc.moc.types import CharNet, PointKind
from rockmoc.thermo.gas_protocol import MOCGas


@dataclass
class AxisymHFResult:
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


@_records_fallbacks
def design_axisym_mln_hf(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    n_chars: int = 40,
    throat_mach: float = 1.01,
    delta: float = 1.0,
    match_exit_area: bool = True,
    store_net: bool = True,
    resample: int | None = None,
) -> AxisymHFResult:
    """Minimum-length nozzle via full MOC net (planar or axisymmetric)."""
    if exit_mach <= 1.0:
        raise DesignError("exit_mach must be > 1")
    if throat_radius <= 0.0:
        raise DesignError("throat_radius must be positive")
    if n_chars < 4:
        raise DesignError("n_chars must be ≥ 4")

    Rt = float(throat_radius)
    Me = float(exit_mach)
    nu_e = float(gas.prandtl_meyer(Me))
    theta_max = 0.5 * nu_e

    def _net(th: float, transition: bool = True) -> CharNet:
        return build_mln_full_net(
            gas, Rt=Rt, theta_max=th, n_chars=n_chars, throat_mach=throat_mach, delta=delta,
            transition=transition,
        )[0]

    if abs(delta) > 0.5:
        # Axisymmetric: sources raise axis ν above 2θ_max, so solve for the
        # corner angle that puts the kernel's last axis point exactly at ν_e.
        from scipy.optimize import brentq

        def _resid(th: float) -> float:
            # probes may wander into over-expanded corners: their fallbacks are not
            # part of the final design, so don't record them
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", MOCFallbackWarning)
                try:
                    return _net(th, transition=False).axis[-1].nu - nu_e
                except (MeshError, ThermoError):
                    return nu_e  # corner too steep: over-expands past ν_max

        # θ_max/ν_e ≈ 0.22–0.24 for γ = 1.2–1.4: try a tight bracket first
        for lo, hi in ((0.20 * nu_e, 0.27 * nu_e), (0.05 * nu_e, 0.5 * nu_e)):
            if _resid(lo) < 0.0 < _resid(hi):
                theta_max = float(brentq(_resid, lo, hi, xtol=1e-10))
                break
        else:
            warn_fallback("θ_max not bracketed; planar value ν_e/2 used")

    net = _net(theta_max)

    # Wall angles should start near θ_max: set throat wall node angle
    if net.wall and len(net.wall) >= 2:
        net.wall[0] = _fp(0.0, Rt, theta_max, theta_max, gas, PointKind.THROAT)

    contour = _wall_to_contour(net.wall, Rt)
    # Force first wall angle to θ_max for contour colouring
    if contour.theta is not None and len(contour.theta):
        th = contour.theta.copy()
        th[0] = theta_max
        th[-1] = 0.0
        contour = Contour(x=contour.x, r=contour.r, theta=th, M=contour.M)

    eps_1d = float(gas.area_ratio(Me))
    if abs(delta) < 0.5:
        r_target = Rt * eps_1d
    else:
        r_target = Rt * math.sqrt(eps_1d)

    if match_exit_area and contour.r[-1] > Rt + 1e-15:
        scale = (r_target - Rt) / (float(contour.r[-1]) - Rt)
        contour = _isotropic_match_exit(contour, Rt, r_target)
        _scale_net(net, Rt, scale, x0=0.0)
        if net.wall:
            net.wall[0] = _fp(0.0, Rt, theta_max, theta_max, gas, PointKind.THROAT)
            last = net.wall[-1]
            net.wall[-1] = _fp(last.x, r_target, 0.0, last.nu, gas, PointKind.EXIT)

    n_out = resample if resample is not None else max(80, n_chars * 3)
    contour_out = contour.resample(n_out)

    if abs(delta) < 0.5:
        eps_g = float(contour_out.r[-1] / Rt)
    else:
        eps_g = float((contour_out.r[-1] / Rt) ** 2)

    axis_M = net.axis[-1].M if net.axis else None
    max_ang = float(np.max(np.abs(contour_out.theta))) if contour_out.theta is not None else theta_max

    return AxisymHFResult(
        kind="mln",
        exit_mach_design=Me,
        gamma=float(gas.gamma),
        throat_radius=Rt,
        n_chars=n_chars,
        theta_max_rad=theta_max,
        contour=contour_out,
        net=net if store_net else None,
        area_ratio_geometric=eps_g,
        area_ratio_1d=eps_1d,
        length=float(contour_out.length),
        max_wall_angle_rad=max_ang,
        axis_exit_mach=float(axis_M) if axis_M is not None else None,
        meta={
            "delta": delta,
            "match_exit_area": match_exit_area,
            "nu_e_rad": nu_e,
            "throat_mach": throat_mach,
            "n_wall_nodes": len(net.wall),
            "n_net_points": len(net.points),
            "algorithm": "full_moc_net_mln",
        },
    )


@_records_fallbacks
def design_axisym_ideal_hf(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    n_chars: int = 40,
    throat_rc_over_rt: float = 1.5,
    throat_mach: float = 1.01,
    delta: float = 1.0,
    match_exit_area: bool = True,
    store_net: bool = True,
) -> AxisymHFResult:
    """Ideal (uniform, axial exit) contour with a circular-arc throat — full MOC.

    The arc attachment angle θ_a is solved so that the last expansion C−
    reaches the axis at M = Me; the downstream wall closes the throat mass flow.
    """
    if exit_mach <= 1.0:
        raise DesignError("exit_mach must be > 1")
    if throat_radius <= 0.0 or throat_rc_over_rt <= 0.0:
        raise DesignError("throat_radius and throat_rc_over_rt must be positive")
    Rt = float(throat_radius)
    Rc = throat_rc_over_rt * Rt
    Me = float(exit_mach)
    nu_e = float(gas.prandtl_meyer(Me))

    net = net_lo = None
    for start in ("hall", "sauer", "uniform"):
        if start == "uniform":
            warn_fallback("Hall/Sauer start not usable (small Rc / low Me); uniform near-sonic start used")
        try:
            kw = dict(gas=gas, Rt=Rt, Rc=Rc, nu_target=nu_e, n_chars=n_chars,
                      throat_mach=throat_mach, delta=delta, initial_line=start)
            net, net_lo = build_ideal_nets(**kw, offsets=(0, 1))  # type: ignore[arg-type]
            break
        except (MeshError, ThermoError):
            net = net_lo = None
    if net is None or net_lo is None:
        raise DesignError(
            f"ideal contour not found for Me={Me}; "
            "try a larger throat_rc_over_rt or lower throat_mach"
        )
    theta_a = _attach_angle(net, Rt, Rc)

    # The attachment is discrete: blend the two bracketing designs linearly in ν_E
    nu_hi, nu_lo = net.axis[-1].nu, net_lo.axis[-1].nu
    wgt = min(max((nu_e - nu_lo) / max(nu_hi - nu_lo, 1e-12), 0.0), 1.0)
    contour = blend_contours(_wall_to_contour(net_lo.wall, Rt), _wall_to_contour(net.wall, Rt),
                             wgt, n=max(200, 6 * n_chars))
    theta_a = (1.0 - wgt) * _attach_angle(net_lo, Rt, Rc) + wgt * theta_a
    if wgt < 0.5:  # keep the bracketing net nearer the blended design
        net = net_lo
    eps_1d = float(gas.area_ratio(Me))
    planar = abs(delta) < 0.5
    r_t = Rt * (eps_1d if planar else math.sqrt(eps_1d))
    if match_exit_area and contour.r[-1] > Rt:
        contour = _isotropic_match_exit(contour, Rt, r_t)
    contour = contour.resample(max(100, n_chars * 3))
    ratio = float(contour.r[-1] / Rt)
    eps_g = ratio if planar else ratio**2
    max_ang = float(np.max(np.abs(contour.theta)))

    return AxisymHFResult(
        kind="ideal",
        exit_mach_design=Me,
        gamma=float(gas.gamma),
        throat_radius=Rt,
        n_chars=n_chars,
        theta_max_rad=theta_a,
        contour=contour,
        net=net if store_net else None,
        area_ratio_geometric=eps_g,
        area_ratio_1d=eps_1d,
        length=float(contour.length),
        max_wall_angle_rad=max_ang,
        axis_exit_mach=float(net.axis[-1].M),
        meta={
            "delta": delta,
            "nu_e_rad": nu_e,
            "throat_mach": throat_mach,
            "throat_rc_over_rt": throat_rc_over_rt,
            "theta_attach_deg": math.degrees(theta_a),
            "initial_line": start,
            "match_exit_area": match_exit_area,
            "n_wall_nodes": len(net.wall),
            "algorithm": "full_moc_net_ideal_arc",
        },
    )


@_records_fallbacks
def design_axisym_tic_hf(
    gas: MOCGas,
    *,
    throat_radius: float,
    exit_mach: float,
    truncate_fraction: float = 0.80,
    n_chars: int = 40,
    **kwargs,
) -> AxisymHFResult:
    """Truncated ideal contour at a fraction of ideal length."""
    if not (0.05 < truncate_fraction <= 1.0):
        raise DesignError("truncate_fraction must be in (0.05, 1]")
    store = kwargs.pop("store_net", True)
    full = design_axisym_ideal_hf(
        gas,
        throat_radius=throat_radius,
        exit_mach=exit_mach,
        n_chars=n_chars,
        store_net=store,
        **kwargs,
    )
    contour = full.contour.truncate(truncate_fraction)
    ratio = float(contour.r[-1] / throat_radius)
    eps_g = ratio if abs(full.meta.get("delta", 1.0)) < 0.5 else ratio**2

    net = full.net
    if net is not None and store:
        x_cut = float(contour.x[-1])
        net.wall = [p for p in net.wall if p.x <= x_cut + 1e-12]
        if net.wall:
            th_e = float(contour.theta[-1]) if contour.theta is not None else net.wall[-1].theta
            Me_e = float(contour.M[-1]) if contour.M is not None else net.wall[-1].M
            nu_e = float(gas.prandtl_meyer(max(Me_e, 1.01)))
            net.wall[-1] = _fp(x_cut, float(contour.r[-1]), th_e, nu_e, gas, PointKind.EXIT)
        net.axis = [p for p in net.axis if p.x <= x_cut + 1e-9]
        net.interior = [[p for p in row if p.x <= x_cut + 1e-9] for row in net.interior]
        net.interior = [row for row in net.interior if len(row) >= 2]
        net.points = [p for p in net.points if p.x <= x_cut + 1e-9]

    max_ang = float(np.max(np.abs(contour.theta))) if contour.theta is not None else full.max_wall_angle_rad
    return AxisymHFResult(
        kind="tic",
        exit_mach_design=full.exit_mach_design,
        gamma=full.gamma,
        throat_radius=throat_radius,
        n_chars=n_chars,
        theta_max_rad=full.theta_max_rad,
        contour=contour,
        net=net,
        area_ratio_geometric=eps_g,
        area_ratio_1d=full.area_ratio_1d,
        length=float(contour.length),
        max_wall_angle_rad=max_ang,
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
    match_exit_area: bool = True,
) -> list[dict]:
    if n_list is None:
        n_list = [12, 20, 40, 60, 80]
    rows: list[dict] = []
    prev = None
    for n_c in n_list:
        r = design_axisym_mln_hf(
            gas,
            throat_radius=throat_radius,
            exit_mach=exit_mach,
            n_chars=n_c,
            delta=delta,
            match_exit_area=match_exit_area,
            store_net=False,
        )
        s = r.summary()
        row = {
            "n_chars": n_c,
            "length_over_Rt": s["length_over_Rt"],
            "area_ratio_geometric": s["area_ratio_geometric"],
            "area_ratio_1d": s["area_ratio_1d"],
            "area_ratio_rel_error": s["area_ratio_rel_error"],
            "theta_max_deg": s["theta_max_deg"],
            "axis_exit_mach": s["axis_exit_mach"],
        }
        if prev is not None:
            c0 = prev.contour.resample(100)
            c1 = r.contour.resample(100)
            r1i = np.interp(c0.x, c1.x, c1.r, left=c1.r[0], right=c1.r[-1])
            row["contour_L2_vs_prev"] = float(
                np.sqrt(np.mean((r1i - c0.r) ** 2)) / throat_radius
            )
        rows.append(row)
        prev = r
    return rows
