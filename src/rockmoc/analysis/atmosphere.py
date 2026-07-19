"""Simple Earth atmosphere pressure model for off-design / aerospike studies.

Valid for approximate engineering use to ~80 km.
"""

from __future__ import annotations

import math

P0 = 101325.0  # Pa
T0 = 288.15  # K
G0 = 9.80665
R_AIR = 287.05
LAPSE = 0.0065  # K/m troposphere


def pressure_at_altitude(h_m: float) -> float:
    """Ambient pressure [Pa] vs geometric altitude [m] AMSL."""
    h = max(float(h_m), 0.0)
    if h <= 11000.0:
        T = T0 - LAPSE * h
        return P0 * (T / T0) ** (G0 / (LAPSE * R_AIR))
    # 11 km state
    T11 = T0 - LAPSE * 11000.0
    p11 = P0 * (T11 / T0) ** (G0 / (LAPSE * R_AIR))
    if h <= 25000.0:
        return p11 * math.exp(-G0 * (h - 11000.0) / (R_AIR * T11))
    # 25 km state
    p25 = p11 * math.exp(-G0 * (25000.0 - 11000.0) / (R_AIR * T11))
    if h <= 80000.0:
        return p25 * math.exp(-G0 * (h - 25000.0) / (R_AIR * 220.0))
    return 0.0


def altitude_from_pressure(pa: float) -> float:
    """Approximate inverse of :func:`pressure_at_altitude` for pa > 0."""
    if pa <= 0:
        return 80_000.0
    if pa >= P0:
        return 0.0
    lo, hi = 0.0, 80_000.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if pressure_at_altitude(mid) > pa:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
