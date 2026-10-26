"""Boundary-layer displacement correction and skin-friction loss for MOC walls.

Primary reference (local copy ``refs/NASA-TM-X-2343.pdf``)
----------------------------------------------------------
Goldman & Vanco, NASA TM-X-2343 (1971): design sharp-edged-throat nozzles by
MOC, then **add displacement thickness** to ideal coordinates. Curvature
neglected.

Method (``method="integral"``, default)
---------------------------------------
Compressible von Kármán momentum-integral equation along the wall arc length
s, with the edge state taken isentropic from the local wall Mach number::

    dθ/ds = Cf/2 − θ [ (2 + H − M²) / (1 + (γ−1)/2 M²) · (1/M) dM/ds  +  δ (1/r) dr/ds ]

(δ = 1 axisymmetric, 0 planar.) Pressure-gradient and history effects enter
through the bracket; the closures are

* skin friction — Eckert/Sommer–Short reference temperature
  ``T*/T_e = 1 + 0.032 M² + 0.58 (T_w/T_e − 1)`` applied to the incompressible
  laws Cf = 0.0256 Re_θ^−1/4 (turbulent, 1/7 power) or Cf = 0.441/Re_θ
  (laminar, Blasius): ``Cf = (ρ*/ρ_e) Cf_inc(Re_θ*)``;
* shape factor — Crocco–Busemann temperature profile in transformed
  coordinates gives exactly ``H = H_i T_w/T_e + r (γ−1)/2 M²``
  (H_i = 9/7 turbulent, 2.59 laminar);
* transport — Bartz (1957) ``μ = 1.184e-7 𝓜^0.5 T^0.6`` [Pa·s, K, kg/kmol]
  and ``Pr = 4γ/(9γ − 5)`` for combustion products.

Fallback (``method="empirical"`` or when the ODE fails)
-------------------------------------------------------
The zero-pressure-gradient closed-form solution of the same closure
(flat-plate correlation at the local edge state). It is what the integral
solver reduces to for dM/ds = dr/ds = 0, so both paths agree there.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from rockmoc.exceptions import DesignError, warn_fallback
from rockmoc.geometry.contour import Contour
from rockmoc.thermo.gas_protocol import MOCGas

R_UNIVERSAL = 8314.462618  # J/(kmol K)

# Incompressible closure Cf_inc = a · Re_θ^(−n) and transformed shape factor H_i
_CLOSURE = {True: (0.0256, 0.25, 9.0 / 7.0), False: (0.441, 1.0, 2.59)}


@dataclass(slots=True)
class BLCorrectionResult:
    """Inviscid + BL-corrected contours, BL distributions and friction loss."""

    inviscid: Contour
    corrected: Contour
    delta_star: NDArray[np.floating]  # m
    theta_mom: NDArray[np.floating]  # momentum thickness m
    Re_x: NDArray[np.floating]
    meta: dict
    cf: NDArray[np.floating] | None = None  # τ_w / (½ ρ_e U_e²)
    tau_w: NDArray[np.floating] | None = None  # Pa
    H: NDArray[np.floating] | None = None  # δ*/θ
    phi: NDArray[np.floating] | None = None  # energy thickness m (two-layer path)
    K_accel: NDArray[np.floating] | None = None  # acceleration parameter
    laminar_weight: NDArray[np.floating] | None = None  # 0 turbulent … 1 relaminarized

    @property
    def delta_cf_friction(self) -> float:
        """Thrust-coefficient loss from wall shear, ΔCf = F_τ / (p_c A_t)."""
        return float(self.meta.get("delta_cf_friction", 0.0))

    def summary(self) -> dict:
        return {
            "delta_star_exit_m": float(self.delta_star[-1]),
            "delta_star_max_m": float(np.max(self.delta_star)),
            "theta_exit_m": float(self.theta_mom[-1]),
            "exit_radius_inviscid_m": float(self.inviscid.r[-1]),
            "exit_radius_corrected_m": float(self.corrected.r[-1]),
            "length_m": float(self.corrected.length),
            **self.meta,
        }


def bartz_viscosity(T: ArrayLike, molar_mass: float) -> NDArray | float:
    """Combustion-gas viscosity [Pa·s] (Bartz 1957), T [K], molar mass [kg/kmol]."""
    return 1.184e-7 * math.sqrt(molar_mass) * np.asarray(T, dtype=np.float64) ** 0.6


def bartz_prandtl(gamma: float) -> float:
    """Eucken-type Prandtl number used by Bartz, Pr = 4γ / (9γ − 5)."""
    return 4.0 * gamma / (9.0 * gamma - 5.0)


def _theta_flat_plate(s, rho_e, U_e, rho_ref, mu_ref, a: float, n: float):
    """ZPG solution of dθ/ds = ½ (ρ*/ρ_e) a (ρ* U θ / μ*)^(−n)."""
    k = 0.5 * a * (rho_ref / rho_e) * (mu_ref / (rho_ref * U_e)) ** n
    return ((n + 1.0) * k * s) ** (1.0 / (n + 1.0))


def two_layer_thicknesses(M: float, tw_over_t0: float, delta_ratio: float, gamma: float,
                          n: int = 400) -> tuple[float, float]:
    """Bartz two-thickness profiles → (H = δ*/θ, φ/θ).

    Velocity u/U = (y/δ)^(1/7); stagnation-temperature defect
    (T0 − T0,local)/(T0 − Tw) = 1 − (y/Δ)^(1/7) over its own thickness Δ
    (``delta_ratio`` = Δ/δ). Δ/δ = 1 reproduces the single-layer value; a
    thick cooled thermal layer lowers H, and H < 0 is physical (Bartz 1965).
    """
    T0Te = 1.0 + 0.5 * (gamma - 1.0) * M * M
    t = np.linspace(0.0, max(1.0, delta_ratio) ** (1.0 / 7.0), n)  # y/δ = t^7 resolves the wall
    eta, w = t**7, 7.0 * t**6
    f = np.minimum(eta, 1.0) ** (1.0 / 7.0)
    gdef = 1.0 - np.minimum(eta / delta_ratio, 1.0) ** (1.0 / 7.0)
    rho = 1.0 / ((1.0 - (1.0 - tw_over_t0) * gdef) * T0Te - (T0Te - 1.0) * f * f)
    trapz = getattr(np, "trapezoid", None) or np.trapz
    th = trapz(rho * f * (1.0 - f) * w, t)
    return float(trapz((1.0 - rho * f) * w, t) / th), float(trapz(rho * f * gdef * w, t) / th)


def _delta_ratio(phi_over_theta: float, M: float, tw: float, gamma: float) -> float:
    """Invert φ/θ(Δ/δ) (monotone increasing) for the thermal/velocity thickness ratio."""
    lo, hi = 1e-3, 40.0
    f = lambda d: two_layer_thicknesses(M, tw, d, gamma, 200)[1] - phi_over_theta  # noqa: E731
    if f(lo) >= 0.0:
        return lo
    if f(hi) <= 0.0:
        return hi
    return float(brentq(f, lo, hi, xtol=1e-4))


def correct_wall_bl(
    contour: Contour,
    gas: MOCGas,
    *,
    throat_radius: float,
    pc: float = 5.0e6,
    Tc: float = 3000.0,
    prandtl: float | None = None,
    turbulent: bool = True,
    axisymmetric: bool = True,
    method: Literal["integral", "empirical"] = "integral",
    wall_temperature: float | ArrayLike | None = None,
    viscosity: Callable[[NDArray], NDArray] | None = None,
    theta0: float | None = None,
    phi0: float | None = None,
    relaminarization: bool = True,
) -> BLCorrectionResult:
    """Displacement-thickness wall correction and skin-friction thrust loss.

    With a prescribed (cooled) ``wall_temperature`` the integral method also
    solves Bartz's energy-integral equation (interaction exponent n = 0, NASA
    TN D-6594)::

        dφ/ds = St (T_aw − T_w)/(T0 − T_w) − φ [d ln(ρU)/ds + δ d ln r/ds + d ln(T0 − T_w)/ds]
        St = (C_f(Re_φ)/2) / (1 − 5 √(C_f/2) [1 − Pr + ln(6/(5Pr + 1))])

    and takes H = δ*/θ from :func:`two_layer_thicknesses` with Δ/δ set by φ/θ,
    so a thick thermal layer (φ/θ up to 5–10 in cooled throats) lowers δ*.

    Relaminarization (``relaminarization=True``): where
    K = (ν0/U²) dU/ds + 0.352 ν0/(U r) dr/ds exceeds 2.88e-6 (NASA TN D-4788)
    the closures blend linearly to laminar between K = 2.5e-6 and 4e-6.

    Parameters
    ----------
    contour:
        Wall (x along axis, r radial) starting at the throat, or at the chamber
        inlet to carry the convergent-section history (subsonic Mach is taken
        from ``contour.M`` or the subsonic area-ratio branch upstream of min r).
    gas:
        Gas model (γ, R); γ is taken at the throat-to-exit mean via ``gas.gamma``.
    pc, Tc:
        Chamber stagnation p [Pa], T [K] for edge-state estimates.
    prandtl:
        Gas Prandtl number; default Bartz ``4γ/(9γ−5)``.
    turbulent:
        Turbulent (1/7 power) or laminar (Blasius) closure.
    method:
        ``"integral"`` momentum-integral ODE (falls back to ``"empirical"``
        with a :class:`MOCFallbackWarning` if it fails), or ``"empirical"``
        flat-plate correlation at the local edge state.
    wall_temperature:
        Wall temperature [K], scalar or per contour point. ``None`` = adiabatic
        wall (T_w = T_aw). Regen-cooled walls (≈ 600–1000 K) give thinner δ*.
    viscosity:
        μ(T) [Pa·s]; default :func:`bartz_viscosity` with 𝓜 = R_u / R.
    theta0:
        Momentum thickness at the first contour point [m]. Default: flat-plate
        value for a run length of one throat radius.
    phi0:
        Energy thickness at the first point [m] (cooled walls). Default θ0
        (Δ = δ); a cooled chamber gives φ0 > θ0.
    """
    if contour.x.size < 3:
        raise DesignError("contour too short for BL correction")
    if pc <= 0 or Tc <= 0 or throat_radius <= 0:
        raise DesignError("pc, Tc and throat_radius must be positive")
    if method not in ("integral", "empirical"):
        raise DesignError(f"unknown BL method {method!r}; use 'integral' or 'empirical'")

    x = np.asarray(contour.x, dtype=np.float64)
    r = np.asarray(contour.r, dtype=np.float64)
    n_pts = len(x)
    g = float(gas.gamma)
    Rg = float(gas.R)
    Pr = bartz_prandtl(g) if prandtl is None else float(prandtl)
    if viscosity is None:
        mol = R_UNIVERSAL / Rg
        viscosity = lambda T: bartz_viscosity(T, mol)  # noqa: E731

    # Local edge Mach from contour.M if present, else from area ratio
    i_t = int(np.argmin(r))  # geometric throat
    if contour.M is not None and len(contour.M) == n_pts:
        M_e = np.maximum(np.asarray(contour.M, dtype=np.float64), 1e-3)
    else:
        M_e = np.ones(n_pts)
        for i in range(n_pts):
            ar = max((r[i] / throat_radius) ** 2 if axisymmetric else (r[i] / throat_radius), 1.0)
            try:
                M_e[i] = float(gas.mach_from_area_ratio(ar, branch="subsonic" if i < i_t else "supersonic"))
            except Exception:
                M_e[i] = 1.0
        M_e = np.maximum(np.maximum.accumulate(M_e), 1e-3)

    # Isentropic edge state
    m2 = 0.5 * (g - 1.0) * M_e**2
    T_e = Tc / (1.0 + m2)
    p_e = pc / (1.0 + m2) ** (g / (g - 1.0))
    rho_e = p_e / (Rg * T_e)
    U_e = M_e * np.sqrt(g * Rg * T_e)
    mu_e = np.asarray(viscosity(T_e), dtype=np.float64)

    ds = np.hypot(np.diff(x), np.diff(r))
    s = np.concatenate([[0.0], np.cumsum(ds)])
    if s[-1] <= 0.0:
        raise DesignError("contour has zero length")
    Re_x = np.maximum(rho_e * U_e * np.maximum(s, 1e-12) / mu_e, 1.0)

    r_rec = Pr ** (1.0 / 3.0) if turbulent else math.sqrt(Pr)
    T_aw = T_e * (1.0 + r_rec * m2)
    if wall_temperature is None:
        T_w = T_aw
    else:
        T_w = np.broadcast_to(np.asarray(wall_temperature, dtype=np.float64), (n_pts,)).copy()
        if np.any(T_w <= 0.0):
            raise DesignError("wall_temperature must be positive")
    T_ref = T_e * (1.0 + 0.032 * M_e**2 + 0.58 * (T_w / T_e - 1.0))
    rho_ref = p_e / (Rg * T_ref)
    mu_ref = np.asarray(viscosity(T_ref), dtype=np.float64)

    a, n_exp, H_i = _CLOSURE[bool(turbulent)]
    a_l, n_l, H_il = _CLOSURE[False]
    T0 = Tc

    # Relaminarization weight w ∈ [0, 1] from the acceleration parameter (TN D-4788)
    nu0 = float(np.asarray(viscosity(np.array([T0]))).ravel()[0]) / (pc / (Rg * T0))
    dUds = np.gradient(U_e, s)
    K = nu0 / U_e**2 * dUds + (0.352 * nu0 / (U_e * np.maximum(r, 1e-300)) * np.gradient(r, s) if axisymmetric else 0.0)
    w_lam = np.clip((K - 2.5e-6) / 1.5e-6, 0.0, 1.0) if (relaminarization and turbulent) else np.zeros(n_pts)
    H_single = ((1.0 - w_lam) * H_i + w_lam * H_il) * T_w / T_e + r_rec * m2

    def cf_of(theta, rho_e_, U_, rho_r, mu_r, w=0.0):
        Re_th = np.maximum(rho_r * U_ * theta / mu_r, 1.0)
        turb = (rho_r / rho_e_) * a * Re_th ** (-n_exp)
        return turb if not np.any(w) else (1.0 - w) * turb + w * (rho_r / rho_e_) * a_l * Re_th ** (-n_l)

    if theta0 is None:
        th0 = float(_theta_flat_plate(throat_radius, rho_e[0], U_e[0], rho_ref[0], mu_ref[0], a, n_exp))
    else:
        th0 = float(theta0)
    if th0 <= 0.0:
        raise DesignError("theta0 must be positive")
    ph0 = th0 if phi0 is None else float(phi0)
    if ph0 <= 0.0:
        raise DesignError("phi0 must be positive")

    # Two-layer energy-integral path: prescribed wall temperature clearly below T0
    two_layer = method == "integral" and wall_temperature is not None and bool(np.all(T_w < 0.98 * T0))
    tw0 = T_w / T0
    crocco_pr = 1.0 - Pr + math.log(6.0 / (5.0 * Pr + 1.0))

    def H_two(th: float, ph: float, M: float, tw: float, w: float) -> float:
        h2 = two_layer_thicknesses(M, tw, _delta_ratio(ph / th, M, tw, g), g, 200)[0]
        return h2 if w <= 0.0 else (1.0 - w) * h2 + w * (H_il * tw * (1.0 + 0.5 * (g - 1.0) * M * M)
                                                        + math.sqrt(Pr) * 0.5 * (g - 1.0) * M * M)

    used, reason = method, None
    theta = phi = None
    H = H_single
    if method == "integral":
        # Pressure-gradient and divergence coefficients on the stations
        dlnM = np.gradient(M_e, s) / M_e
        dlnr = np.gradient(r, s) / np.maximum(r, 1e-300) if axisymmetric else np.zeros(n_pts)
        dln_rhoU = np.gradient(rho_e * U_e, s) / (rho_e * U_e)
        dln_dT = -np.gradient(T_w, s) / np.maximum(T0 - T_w, 1e-9)
        drive = (T_aw - T_w) / np.maximum(T0 - T_w, 1e-9)
        stations = np.vstack([rho_e, U_e, rho_ref, mu_ref, H_single, M_e, dlnM, dlnr, w_lam,
                              tw0, dln_rhoU, dln_dT, drive])

        def rhs(si: float, y: NDArray) -> NDArray:
            ro, U_, rr, mr, Hs, M, dM, dr_, w, tw, dlg, dlt, drv = (np.interp(si, s, row) for row in stations)
            th = max(float(y[0]), 1e-300)
            if two_layer:
                ph = max(float(y[1]), 1e-300)
                Hh = H_two(th, ph, M, tw, w)
            else:
                Hh = Hs
            grad = (2.0 + Hh - M * M) / (1.0 + 0.5 * (g - 1.0) * M * M) * dM + dr_
            dth = 0.5 * cf_of(th, ro, U_, rr, mr, w) - th * grad
            if not two_layer:
                return np.array([dth])
            c2 = 0.5 * cf_of(ph, ro, U_, rr, mr, w)
            St = c2 / max(1.0 - 5.0 * math.sqrt(c2) * crocco_pr, 0.2)
            return np.array([dth, St * drv - ph * (dlg + dr_ + dlt)])

        y0 = [th0, ph0] if two_layer else [th0]
        try:
            sol = solve_ivp(rhs, (0.0, float(s[-1])), y0, method="LSODA", t_eval=s,
                            rtol=1e-6, atol=1e-12 * throat_radius, max_step=float(s[-1]) / 50.0)
            if not sol.success:
                reason = sol.message
            elif sol.y.shape[1] != n_pts or not np.all(np.isfinite(sol.y)) or np.any(sol.y <= 0.0):
                reason = "non-physical momentum/energy thickness"
            else:
                theta = sol.y[0]
                if two_layer:
                    phi = sol.y[1]
                    H = np.array([H_two(float(t_), float(p_), float(m_), float(tw_), float(w_))
                                  for t_, p_, m_, tw_, w_ in zip(theta, phi, M_e, tw0, w_lam)])
        except (ValueError, ArithmeticError) as exc:
            reason = str(exc)
        if theta is None:
            used = "empirical"
            warn_fallback(f"BL momentum-integral solve failed ({reason}); flat-plate correlation used")

    if theta is None:
        # Run length s + s0 so the empirical path also honours θ(0) = θ0
        s0 = (th0 / _theta_flat_plate(1.0, rho_e[0], U_e[0], rho_ref[0], mu_ref[0], a, n_exp)) ** (n_exp + 1.0)
        theta = _theta_flat_plate(s + s0, rho_e, U_e, rho_ref, mu_ref, a, n_exp)
        H = H_single

    delta_star = H * theta
    cf = cf_of(theta, rho_e, U_e, rho_ref, mu_ref, w_lam)
    tau_w = 0.5 * cf * rho_e * U_e**2

    if contour.theta is not None:
        alpha = np.asarray(contour.theta, dtype=np.float64)
    else:
        alpha = np.arctan2(np.gradient(r), np.gradient(x))

    # Wall shear acts along the wall tangent; its axial part opposes thrust.
    # Divergent wall only: convergent-section friction belongs to c* efficiency.
    trapz = getattr(np, "trapezoid", None) or np.trapz
    sd = slice(i_t, None)
    if axisymmetric:
        F_tau = float(trapz((tau_w * np.cos(alpha) * 2.0 * math.pi * r)[sd], s[sd]))
        A_t = math.pi * throat_radius**2
    else:  # per unit depth, both walls of a symmetric 2-D nozzle
        F_tau = float(trapz((tau_w * np.cos(alpha) * 2.0)[sd], s[sd]))
        A_t = 2.0 * throat_radius
    d_cf = F_tau / (pc * A_t)

    # Grow wall along local outward normal (-sin α, cos α)
    r_corr = r + delta_star * np.cos(alpha)
    x_corr = x - delta_star * np.sin(alpha)
    for i in range(1, n_pts):
        if x_corr[i] <= x_corr[i - 1]:
            x_corr[i] = x_corr[i - 1] + 1e-12
    r_corr = np.maximum(r_corr, throat_radius * 0.99)

    corrected = Contour(x=x_corr, r=r_corr, theta=alpha, M=M_e)

    return BLCorrectionResult(
        inviscid=contour,
        corrected=corrected,
        delta_star=delta_star,
        theta_mom=theta,
        Re_x=Re_x,
        cf=cf,
        tau_w=tau_w,
        H=H,
        phi=phi,
        K_accel=K,
        laminar_weight=w_lam,
        meta={
            "method": ("momentum_energy_integral" if phi is not None else "momentum_integral")
            if used == "integral" else "empirical_flat_plate",
            "relaminarized": bool(np.any(w_lam > 0.5)),
            "method_requested": method,
            "fallback_reason": reason,
            "turbulent": turbulent,
            "pc_Pa": pc,
            "Tc_K": Tc,
            "prandtl": Pr,
            "wall": "adiabatic" if wall_temperature is None else "prescribed",
            "axisymmetric": axisymmetric,
            "friction_force_N": F_tau,
            "delta_cf_friction": d_cf,
            "reference": "von Kármán momentum integral, Eckert ref. T, Bartz transport; TM-X-2343 δ* rule",
        },
    )
