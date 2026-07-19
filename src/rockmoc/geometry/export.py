"""Additional geometry export helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from rockmoc.geometry.contour import Contour


def export_numpy(contour: Contour, path: str | Path) -> None:
    """Save contour arrays to a ``.npz`` archive."""
    path = Path(path)
    payload = {"x": contour.x, "r": contour.r, "theta": contour.theta}
    if contour.M is not None:
        payload["M"] = contour.M
    np.savez(path, **payload)


def export_dxf_polyline(contour: Contour, path: str | Path) -> None:
    """Write a minimal ASCII DXF LW polyline of the wall (no external deps)."""
    path = Path(path)
    lines: list[str] = [
        "0",
        "SECTION",
        "2",
        "HEADER",
        "0",
        "ENDSEC",
        "0",
        "SECTION",
        "2",
        "ENTITIES",
    ]
    for x, r in zip(contour.x, contour.r, strict=True):
        lines.extend(
            [
                "0",
                "POINT",
                "8",
                "NOZZLE",
                "10",
                f"{x:.9g}",
                "20",
                f"{r:.9g}",
                "30",
                "0.0",
            ]
        )
    # Also emit a LINE strip between consecutive points
    for i in range(len(contour.x) - 1):
        lines.extend(
            [
                "0",
                "LINE",
                "8",
                "NOZZLE",
                "10",
                f"{contour.x[i]:.9g}",
                "20",
                f"{contour.r[i]:.9g}",
                "30",
                "0.0",
                "11",
                f"{contour.x[i + 1]:.9g}",
                "21",
                f"{contour.r[i + 1]:.9g}",
                "31",
                "0.0",
            ]
        )
    lines.extend(["0", "ENDSEC", "0", "EOF"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
