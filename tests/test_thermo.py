"""Tests for perfect-gas isentropic and Prandtl–Meyer relations."""

from __future__ import annotations

import math

import numpy as np
import pytest

from rockmoc.exceptions import ThermoError
from rockmoc.thermo import PerfectGas, mach_from_prandtl_meyer, prandtl_meyer


def test_prandtl_meyer_known_value():
    # ν(M=2, γ=1.4) ≈ 26.3798 deg
    gas = PerfectGas(gamma=1.4)
    nu = float(gas.prandtl_meyer(2.0))
    assert abs(math.degrees(nu) - 26.37976) < 0.01


def test_prandtl_meyer_inverse_roundtrip():
    gas = PerfectGas(gamma=1.4)
    for M in [1.01, 1.5, 2.0, 3.0, 5.0, 8.0]:
        nu = float(gas.prandtl_meyer(M))
        M2 = float(gas.mach_from_prandtl_meyer(nu))
        assert abs(M2 - M) < 1e-8


def test_area_mach_known_values():
    gas = PerfectGas(gamma=1.4)
    # Standard isentropic table values
    assert abs(float(gas.area_ratio(2.0)) - 1.6875) < 1e-3
    assert abs(float(gas.area_ratio(3.0)) - 4.2346) < 1e-2


def test_area_mach_inverse_roundtrip():
    gas = PerfectGas(gamma=1.25)
    for M in [1.2, 2.0, 4.0, 6.0]:
        eps = float(gas.area_ratio(M))
        M2 = float(gas.mach_from_area_ratio(eps, branch="supersonic"))
        assert abs(M2 - M) < 1e-7


def test_subsonic_branch():
    gas = PerfectGas(gamma=1.4)
    eps = float(gas.area_ratio(0.5))
    M = float(gas.mach_from_area_ratio(eps, branch="subsonic"))
    assert abs(M - 0.5) < 1e-6


def test_pressure_temperature_ratios():
    gas = PerfectGas(gamma=1.4)
    M = 2.0
    # p0/p = (1 + 0.2 M^2)^(3.5)
    p0_p = float(gas.pressure_ratio(M))
    assert abs(p0_p - (1.0 + 0.2 * 4) ** 3.5) < 1e-10
    T0_T = float(gas.temperature_ratio(M))
    assert abs(T0_T - (1.0 + 0.2 * 4)) < 1e-12


def test_invalid_gamma():
    with pytest.raises(ThermoError):
        PerfectGas(gamma=0.9)


def test_module_level_helpers():
    nu = float(prandtl_meyer(2.0, 1.4))
    M = float(mach_from_prandtl_meyer(nu, 1.4))
    assert abs(M - 2.0) < 1e-8


def test_from_molecular_weight():
    gas = PerfectGas.from_molecular_weight(1.25, 22.0)
    assert gas.R == pytest.approx(8314.462618 / 22.0)


def test_thrust_coefficient_positive():
    gas = PerfectGas(gamma=1.2)
    Me = 3.0
    pe_pc = float(gas.pe_pc_from_mach(Me))
    Cf = gas.thrust_coefficient(pe_pc, 0.0, exit_mach=Me)
    assert Cf > 1.0


def test_vectorized_pm():
    gas = PerfectGas(gamma=1.4)
    M = np.array([1.5, 2.0, 2.5])
    nu = gas.prandtl_meyer(M)
    assert nu.shape == (3,)
    assert np.all(nu > 0)
