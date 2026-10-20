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


class MOCFallbackWarning(UserWarning):
    """A numerical fallback was used; the result is valid but degraded.

    Designers record these in ``result.meta["warnings"]`` (see
    :func:`collect_fallbacks`) and re-emit them.
    """


def warn_fallback(message: str) -> None:
    import warnings

    warnings.warn(message, MOCFallbackWarning, stacklevel=3)


class collect_fallbacks:  # noqa: N801 - used as a context manager
    """Context manager: capture :class:`MOCFallbackWarning` messages into ``self.messages``.

    Duplicates are collapsed; on exit each distinct message is re-emitted once.
    """

    def __enter__(self) -> list[str]:
        import warnings

        self._cm = warnings.catch_warnings(record=True)
        self._caught = self._cm.__enter__()
        warnings.simplefilter("always", MOCFallbackWarning)
        self.messages: list[str] = []
        return self.messages

    def __exit__(self, *exc: object) -> None:
        import warnings

        self._cm.__exit__(*exc)
        for w in self._caught:
            if issubclass(w.category, MOCFallbackWarning):
                msg = str(w.message)
                if msg not in self.messages:
                    self.messages.append(msg)
            else:
                warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)
        for msg in self.messages:
            warnings.warn(msg, MOCFallbackWarning, stacklevel=2)
