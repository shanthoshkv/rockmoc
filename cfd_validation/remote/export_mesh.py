"""Runs on the PC (~/venv). Saves the 2-D face polygons of each case's mesh -> mesh.npz.

Uses the 'back' patch (wedge or empty side face), so every cell appears once as
a polygon in (x, r).  Usage: python export_mesh.py <case> [<case> ...]
"""
import sys
from pathlib import Path

import numpy as np
import pyvista as pv

for c in sys.argv[1:]:
    case = Path(c)
    try:
        (case / "case.foam").touch()
        rd = pv.OpenFOAMReader(str(case / "case.foam"))
        rd.set_active_time_value(rd.time_values[0])
        rd.disable_all_cell_arrays()
        rd.disable_all_point_arrays()
        rd.enable_all_patch_arrays()
        mb = rd.read()
        back = mb["boundary"]["back"].extract_surface()
        P = np.asarray(back.points)
        r = np.hypot(P[:, 1], P[:, 2]) * np.sign(P[:, 1] + 1e-30)
        f = np.asarray(back.faces)
        np.savez_compressed(case / "mesh.npz", x=P[:, 0].astype(np.float64), r=r, faces=f.astype(np.int64),
                            ncell=mb["internalMesh"].n_cells)
        print(case.name, "ok", back.n_cells)
    except Exception as e:  # noqa: BLE001
        print(case.name, "FAIL", e)
