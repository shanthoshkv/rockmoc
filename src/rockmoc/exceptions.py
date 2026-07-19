"""rockmoc exception hierarchy."""

from __future__ import annotations


class RockmocError(Exception):
    """Base exception for all rockmoc errors."""


class DesignError(RockmocError):
    """Invalid design inputs or failed design construction."""


class ConvergenceError(RockmocError):
    """Numerical solver failed to converge."""


class ThermoError(RockmocError):
    """Invalid thermodynamic state or property evaluation."""


class MeshError(RockmocError):
    """Characteristic mesh construction failure."""
