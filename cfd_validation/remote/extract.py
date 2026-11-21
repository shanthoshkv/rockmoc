"""Runs on the PC (in ~/venv). Latest time of a reconstructed case -> cfd.npz.

Saves internal cells, wall faces, an exit-plane cut, inlet faces and the
inlet/outlet mass-flow history, all as plain arrays.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

case = Path(sys.argv[1])
m = json.loads((case / "meta.json").read_text())
(case / "case.foam").touch()
rd = pv.OpenFOAMReader(str(case / "case.foam"))
rd.set_active_time_value(rd.time_values[-1])
rd.cell_to_point_creation = False
rd.enable_all_patch_arrays()
mb = rd.read()
im = mb["internalMesh"]
bd = mb["boundary"]
planar = m["planar"]
out = {"time": rd.time_values[-1]}


def rad(pts):
    return np.abs(pts[:, 1]) if planar else np.hypot(pts[:, 1], pts[:, 2])


def dump(prefix, mesh, area=False):
    c = mesh.cell_centers().points
    out[prefix + "x"] = c[:, 0]
    out[prefix + "r"] = rad(c)
    for k in mesh.cell_data.keys():
        a = np.asarray(mesh.cell_data[k])
        if a.ndim == 2 and a.shape[1] == 3:
            out[prefix + k + "_x"] = a[:, 0]
            rr = rad(c)
            if planar:
                out[prefix + k + "_r"] = a[:, 1] * np.sign(c[:, 1])
            else:
                out[prefix + k + "_r"] = (a[:, 1] * c[:, 1] + a[:, 2] * c[:, 2]) / np.maximum(rr, 1e-30)
        elif a.ndim == 1:
            out[prefix + k] = a
    if area:
        s = mesh.compute_cell_sizes(length=False, volume=False)
        out[prefix + "A"] = np.asarray(s.cell_data["Area"])
        n = mesh.compute_normals(cell_normals=True, point_normals=False) if hasattr(mesh, "compute_normals") else None
        if n is not None:
            nn = np.asarray(n.cell_data["Normals"])
            out[prefix + "nx"] = nn[:, 0]


dump("c_", im)
dump("w_", bd["wall"].extract_surface(), area=True)
dump("i_", bd["inlet"].extract_surface(), area=True)
xe = m["wall_x"][-1]
L = m["wall_x"][-1] - m["wall_x"][0]
cut = im.slice(normal="x", origin=(xe - 1e-4 * L, 0, 0))
dump("e_", cut, area=True)

for name in ("flowRates", "flowRateOut"):
    fs = sorted((case / "postProcessing" / name).glob("*/surfaceFieldValue.dat"))
    if fs:
        rows = [l.split() for f in fs for l in f.read_text().splitlines() if l and not l.startswith("#")]
        out["h_" + name] = np.array(rows, dtype=float)

np.savez_compressed(case / "cfd.npz", **out)
print("extract ok", {k: np.shape(v) for k, v in out.items() if k.startswith(("c_x", "w_x", "e_x"))})
