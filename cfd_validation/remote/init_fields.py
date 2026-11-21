"""Runs on the PC (in ~/venv). Writes a 1-D isentropic initial field into 0/p,T,U.

Reads meta.json and 0/C (from `postProcess -func writeCellCentres`). Each 0/<f> file
carries a line `internalField uniform ...;` that is replaced by a nonuniform list.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

case = Path(sys.argv[1])
m = json.loads((case / "meta.json").read_text())
g = m["gamma_init"]


def read_vectors(text, start):
    k = text.index("List<vector>", start)
    i = text.index("(", k)
    n = int(text[k + len("List<vector>"):i].split()[-1])
    body = text[i + 1:]
    pts = re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", body[: body.index("\n)")])
    a = np.array(pts, dtype=float)
    assert n is None or len(a) == n, (n, len(a))
    return a


Ctext = (case / "0" / "C").read_text()
C = read_vectors(Ctext, Ctext.index("internalField"))
x = C[:, 0]
r = np.abs(C[:, 1]) if m["planar"] else np.hypot(C[:, 1], C[:, 2])

xw, rw = np.array(m["wall_x"]), np.array(m["wall_r"])
Rt = m["Rt"]


def area_ratio(M):
    return (1 / M) * ((2 / (g + 1)) * (1 + 0.5 * (g - 1) * M * M)) ** ((g + 1) / (2 * (g - 1)))


def mach_from_eps(eps, sup):
    lo = np.where(sup, 1.0, 1e-4)
    hi = np.where(sup, 50.0, 1.0)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f = area_ratio(mid) - eps
        # subsonic: A/A* decreases with M; supersonic: increases
        up = np.where(sup, f < 0, f > 0)
        lo = np.where(up, mid, lo)
        hi = np.where(up, hi, mid)
    return 0.5 * (lo + hi)


if m["mode"] == "throat":
    Mi, pi, Ti = m["inlet_M"], m["inlet_p"], m["inlet_T"]
    f = 1 + 0.5 * (g - 1) * Mi * Mi
    p0, T0 = pi * f ** (g / (g - 1)), Ti * f
else:
    p0, T0 = m["p0"], m["T0"]

xe = xw[-1]
xc = np.clip(x, xw[0], xe)
rwall = np.interp(xc, xw, rw)
eps = (rwall / Rt) ** (1 if m["planar"] else 2)
eps = np.maximum(eps, 1.0 + 1e-9)
M = mach_from_eps(eps, x >= 0.0)
if m["mode"] == "throat":
    M = np.maximum(M, m["inlet_M"])
tt = 1 + 0.5 * (g - 1) * M * M
T = T0 / tt
p = p0 * tt ** (-g / (g - 1))
U = M * np.sqrt(g * m["R_init"] * T)
if m.get("plume"):
    pl = m["plume"]
    out = x > xe
    Me = M[out]
    core = out & (r <= rw[-1])
    amb = out & (r > rw[-1])
    # core keeps exit state (computed above with x clipped to exit); ambient at rest
    p[amb], T[amb], U[amb] = pl["pa"], pl["Ta"], 0.0
    p[core] = np.maximum(p[core], 0.3 * pl["pa"])
    _ = Me

Ux = U.copy()
fields = {
    "p": "\n".join(f"{v:.8g}" for v in p),
    "T": "\n".join(f"{v:.8g}" for v in T),
    "U": "\n".join(f"({v:.8g} 0 0)" for v in Ux),
}
types = {"p": "scalar", "T": "scalar", "U": "vector"}
N = len(x)
for name, body in fields.items():
    f = case / "0" / name
    s = f.read_text()
    s = re.sub(r"internalField\s+uniform[^;]*;",
               lambda _: f"internalField nonuniform List<{types[name]}>\n{N}\n(\n{body}\n)\n;", s, count=1)
    f.write_text(s)

# prescribed wall temperature profile (x [m] -> Tw [K])
if m.get("Tw_profile"):
    xs, Tw = np.array(m["Tw_profile"][0]), np.array(m["Tw_profile"][1])
    j = Ctext.index("wall", Ctext.index("boundaryField"))
    Cw = read_vectors(Ctext, j)
    vals = np.interp(Cw[:, 0], xs, Tw)
    f = case / "0" / "T"
    s = f.read_text()
    lst = "\n".join(f"{v:.6g}" for v in vals)
    s = s.replace("value uniform TWALL;", f"value nonuniform List<scalar>\n{len(vals)}\n(\n{lst}\n)\n;")
    f.write_text(s)

print(f"init: {N} cells, p0={p0:.6g} T0={T0:.6g}, M range {M.min():.3f}-{M.max():.3f}")
