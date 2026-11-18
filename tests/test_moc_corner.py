"""Axisym MLN corner: wall angle past the corner exceeds θ_max physically, not numerically."""

from __future__ import annotations

import math

import numpy as np

from rockmoc.moc._common import _fp
from rockmoc.moc.axisym_hf import design_axisym_mln_hf
from rockmoc.moc.mesh import interior_point
from rockmoc.moc.types import PointKind
from rockmoc.thermo.perfect_gas import PerfectGas


def test_axisym_unit_process_matches_exact_source_flow():
    """Second-order convergence of the axisym interior point to exact conical source flow."""
    g = PerfectGas(1.4)
    ar0 = float(g.area_ratio(2.0))

    def exact(x, r):
        M = float(g.mach_from_area_ratio(ar0 * (x * x + r * r), branch="supersonic"))
        return math.atan2(r, x), float(g.prandtl_meyer(M))

    errs = []
    for n in (20, 40):
        phi = np.radians(5.0 + 20.0 * np.arange(n + 1) / n)
        line = [_fp(math.cos(p), math.sin(p), *exact(math.cos(p), math.sin(p)), g,
                    PointKind.INTERIOR) for p in phi]
        err = 0.0
        for _ in range(n // 2):
            line = [interior_point(line[i + 1], line[i], g) for i in range(len(line) - 1)]
            err = max(err, *(max(abs(p.theta - exact(p.x, p.r)[0]),
                                 abs(p.nu - exact(p.x, p.r)[1])) for p in line))
        errs.append(err)
    assert errs[0] < math.radians(2e-3)
    assert errs[0] / errs[1] > 3.5  # ~4 = second order


def test_axisym_mln_wall_peak_past_corner_is_converged():
    """γ=1.4, Me=3: corner θ_max ≈ 11.7° but wall peaks ≈ 15.5° downstream (axisym relief).

    The peak is mesh-converged and the mass-streamline slope equals the
    characteristic θ, so it is the physical axisym MLN wall, not overshoot.
    """
    gas = PerfectGas(1.4)
    nu_e = float(gas.prandtl_meyer(3.0))
    peaks = []
    for n in (16, 32, 64):
        r = design_axisym_mln_hf(gas, throat_radius=1.0, exit_mach=3.0, n_chars=n,
                                 match_exit_area=False)
        assert abs(math.degrees(r.theta_max_rad) - 11.73) < 0.05
        assert r.theta_max_rad < 0.5 * nu_e
        w = r.net.wall  # raw MOC wall nodes (skip forced throat node)
        x, rr, th = (np.array([getattr(p, f) for p in w[1:]]) for f in ("x", "r", "theta"))
        th_g = np.degrees(np.arctan2(np.diff(rr), np.diff(x)))
        th_s = np.degrees(0.5 * (th[1:] + th[:-1]))
        assert np.max(np.abs(th_g[:-1] - th_s[:-1])) < 4.0 / n  # O(1/n) at row seam
        peaks.append(float(np.max(np.degrees(th))))
        assert abs(w[-2].theta) < math.radians(4.0) / n  # exit kink → 0 with n
    assert peaks[0] > 15.0
    assert abs(peaks[2] - peaks[1]) < abs(peaks[1] - peaks[0]) + 0.05
    assert abs(peaks[2] - peaks[1]) < 0.1
