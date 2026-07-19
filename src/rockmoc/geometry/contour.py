"""Wall contour representation and export."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from rockmoc.exceptions import DesignError


@dataclass(slots=True)
class Contour:
    """Axisymmetric wall contour polyline.

    ``x`` and ``r`` are SI meters. ``theta`` is wall geometric angle [rad]
    (optional; estimated from finite differences if not provided).
    """

    x: NDArray[np.floating]
    r: NDArray[np.floating]
    theta: NDArray[np.floating] | None = None
    M: NDArray[np.floating] | None = None

    def __post_init__(self) -> None:
        self.x = np.asarray(self.x, dtype=np.float64).ravel()
        self.r = np.asarray(self.r, dtype=np.float64).ravel()
        if self.x.size != self.r.size:
            raise DesignError("Contour x and r must have the same length")
        if self.x.size < 2:
            raise DesignError("Contour needs at least 2 points")
        if self.theta is not None:
            self.theta = np.asarray(self.theta, dtype=np.float64).ravel()
        else:
            self.theta = self._estimate_theta()
        if self.M is not None:
            self.M = np.asarray(self.M, dtype=np.float64).ravel()

    def _estimate_theta(self) -> NDArray[np.floating]:
        dx = np.gradient(self.x)
        dr = np.gradient(self.r)
        return np.arctan2(dr, dx)

    @property
    def length(self) -> float:
        return float(self.x[-1] - self.x[0])

    @property
    def exit_radius(self) -> float:
        return float(self.r[-1])

    @property
    def throat_radius(self) -> float:
        return float(np.min(self.r))

    def area_ratio(self, throat_radius: float | None = None) -> float:
        Rt = self.throat_radius if throat_radius is None else throat_radius
        return float((self.exit_radius / Rt) ** 2)

    def resample(self, n: int) -> Contour:
        """Arc-length resample to ``n`` points."""
        if n < 2:
            raise DesignError("n must be ≥ 2")
        dx = np.diff(self.x)
        dr = np.diff(self.r)
        seg = np.hypot(dx, dr)
        s = np.concatenate([[0.0], np.cumsum(seg)])
        if s[-1] <= 0.0:
            raise DesignError("Contour has zero length")
        s_new = np.linspace(0.0, s[-1], n)
        x = np.interp(s_new, s, self.x)
        r = np.interp(s_new, s, self.r)
        th = np.interp(s_new, s, self.theta) if self.theta is not None else None
        M = np.interp(s_new, s, self.M) if self.M is not None else None
        return Contour(x=x, r=r, theta=th, M=M)

    def truncate(self, fraction: float) -> Contour:
        """Keep the first ``fraction`` of axial length (0–1]."""
        if not (0.0 < fraction <= 1.0):
            raise DesignError("truncate fraction must be in (0, 1]")
        if fraction == 1.0:
            return Contour(x=self.x.copy(), r=self.r.copy(), theta=self.theta.copy() if self.theta is not None else None, M=self.M.copy() if self.M is not None else None)
        x_end = self.x[0] + fraction * (self.x[-1] - self.x[0])
        # include points up to x_end plus interpolated endpoint
        mask = self.x <= x_end + 1e-15
        if not np.any(mask):
            raise DesignError("truncate removed all points")
        x = self.x[mask].copy()
        r = self.r[mask].copy()
        th = self.theta[mask].copy() if self.theta is not None else None
        M = self.M[mask].copy() if self.M is not None else None
        if x[-1] < x_end - 1e-15:
            # interpolate one more point
            x_next = self.x[len(x)]
            r_next = self.r[len(x)]
            t = (x_end - x[-1]) / (x_next - x[-1])
            x = np.append(x, x_end)
            r = np.append(r, r[-1] + t * (r_next - r[-1]))
            if th is not None:
                th_next = self.theta[len(th)]
                th = np.append(th, th[-1] + t * (th_next - th[-1]))
            if M is not None and self.M is not None:
                M_next = self.M[len(M)]
                M = np.append(M, M[-1] + t * (M_next - M[-1]))
        return Contour(x=x, r=r, theta=th, M=M)

    def export_csv(self, path: str | Path) -> None:
        path = Path(path)
        header = "x_m,r_m,theta_rad"
        data = [self.x, self.r, self.theta]
        if self.M is not None:
            header += ",M"
            data.append(self.M)
        arr = np.column_stack(data)
        np.savetxt(path, arr, delimiter=",", header=header, comments="")

    def export_json(self, path: str | Path, extra: dict[str, Any] | None = None) -> None:
        path = Path(path)
        payload: dict[str, Any] = {
            "x_m": self.x.tolist(),
            "r_m": self.r.tolist(),
            "theta_rad": self.theta.tolist() if self.theta is not None else None,
            "M": self.M.tolist() if self.M is not None else None,
        }
        if extra:
            payload["meta"] = extra
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    @classmethod
    def from_wall_points(cls, wall: list[Any]) -> Contour:
        x = np.array([p.x for p in wall], dtype=np.float64)
        r = np.array([p.r for p in wall], dtype=np.float64)
        th = np.array([p.theta for p in wall], dtype=np.float64)
        M = np.array([p.M for p in wall], dtype=np.float64)
        # Ensure strictly increasing x for export/resample (stable sort)
        order = np.argsort(x, kind="mergesort")
        x, r, th, M = x[order], r[order], th[order], M[order]
        # Drop near-duplicate x
        keep = np.ones(len(x), dtype=bool)
        keep[1:] = np.diff(x) > 1e-15 * (1.0 + np.abs(x[1:]))
        return cls(x=x[keep], r=r[keep], theta=th[keep], M=M[keep])
