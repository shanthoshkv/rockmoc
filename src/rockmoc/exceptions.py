"""rockmoc exception hierarchy with student-friendly tips."""

from __future__ import annotations


class RockmocError(Exception):
    """Base exception for all rockmoc errors."""

    def __init__(
        self,
        message: str,
        *,
        tip: str | None = None,
        example: str | None = None,
    ) -> None:
        self.user_message = message
        self.tip = tip
        self.example = example
        parts = [message]
        if tip:
            parts.append(f"  Fix: {tip}")
        if example:
            parts.append(f"  Example: {example}")
        super().__init__("\n".join(parts))


class DesignError(RockmocError):
    """Invalid design inputs or failed design construction."""


class ConvergenceError(RockmocError):
    """Numerical solver failed to converge."""


class ThermoError(RockmocError):
    """Invalid thermodynamic state or property evaluation."""


class MeshError(RockmocError):
    """Characteristic mesh construction failure."""


def design_error(
    message: str,
    *,
    tip: str | None = None,
    example: str | None = None,
) -> DesignError:
    """Build a DesignError with optional Fix/Example lines for students."""
    return DesignError(message, tip=tip, example=example)
