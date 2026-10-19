"""Unit tests for MOC point solvers (planar invariants)."""

from __future__ import annotations

import math

from rockmoc.moc.mesh import axis_point, interior_point
from rockmoc.moc.types import FlowPoint, PointKind
from rockmoc.thermo import PerfectGas


def test_planar_interior_invariants():
    """For δ=0, θ±ν are conserved along C±."""
    gas = PerfectGas(gamma=1.4)
    # Point A (upper) on C−: θ=10°, M=2.0; point B (lower) on C+: θ=5°, M=1.8
    a = FlowPoint.from_mach_theta(
        x=0.0, r=0.8, M=2.0, theta=math.radians(10.0), gas=gas, kind=PointKind.INTERIOR
    )
    b = FlowPoint.from_mach_theta(
        x=0.1, r=0.5, M=1.8, theta=math.radians(5.0), gas=gas, kind=PointKind.INTERIOR
    )
    p = interior_point(a, b, gas, delta=0.0)
    # Planar: (θ+ν) conserved along C− from A, (θ−ν) along C+ from B
    assert abs((p.theta + p.nu) - (a.theta + a.nu)) < 1e-8
    assert abs((p.theta - p.nu) - (b.theta - b.nu)) < 1e-8
    assert p.x > max(a.x, b.x)
    assert p.r > 0
    assert p.M >= 1.0


def test_axis_point_theta_zero():
    gas = PerfectGas(gamma=1.4)
    b = FlowPoint.from_mach_theta(
        x=0.2, r=0.3, M=2.0, theta=math.radians(8.0), gas=gas
    )
    p = axis_point(b, gas, delta=0.0)
    assert p.r == 0.0
    assert abs(p.theta) < 1e-12
    assert p.M >= 1.0
    assert p.x > b.x  # downstream along C−
    assert abs(p.nu - (b.theta + b.nu)) < 1e-8


def test_axisymmetric_interior_runs():
    gas = PerfectGas(gamma=1.25)
    a = FlowPoint.from_mach_theta(0.0, 0.4, 1.6, math.radians(4.0), gas)
    b = FlowPoint.from_mach_theta(0.05, 0.6, 1.9, math.radians(7.0), gas)
    p = interior_point(a, b, gas, delta=1.0)
    assert p.M > 1.0
    assert p.r >= 0.0
    assert math.isfinite(p.x)
