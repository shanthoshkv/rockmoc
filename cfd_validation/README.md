# CFD validation of rockmoc (OpenFOAM)

Automated rockmoc-vs-CFD comparison. Each case designs a nozzle with rockmoc and meshes the wall as a 2.5° wedge (planar slab for 2-D). It then runs OpenFOAM v2412 `rhoCentralFoamAxi` on the remote PC and compares the result against what rockmoc predicts.

```
py -3.12 cfd_validation/cfdval.py build   [case ...]   # write cases/<name> (OpenFOAM case + ref.json)
py -3.12 cfd_validation/cfdval.py launch  [case ...]   # push to PC, 2 lanes x 8 MPI ranks, detached
py -3.12 cfd_validation/cfdval.py status
py -3.12 cfd_validation/cfdval.py pull                 # cfd.npz + logs of finished cases -> results/raw
py -3.12 cfd_validation/cfdval.py compare              # -> results/REPORT.md, results/figures/*.png
```

Remote side (`remote/`): `init_fields.py` writes a 1-D isentropic initial field. `extract.py` (pyvista) reduces the result to cells, wall faces, an exit-plane cut, inlet faces and the mass-flow history. Each case's `Allrun` is restartable: `DONE` marks a case finished, `FAILED` names the step that failed.

| tier | what is compared |
|---|---|
| T1 | inviscid, perfect gas, supersonic inlet on the MOC start line: MOC net field (M, θ), axis M/p/T/ρ/u, wall M, exit uniformity, Isp. Covers every designer: MLN (Me, γ, n_char sweeps), ideal, TIC, Rao TOP/TOC/control-surface, cone, production, planar |
| GCI | coarse/medium/fine (√2 per direction) on MLN, Rao 80 %, cone 15° |
| T2 | full nozzles with chamber + converging section: Cd and throat-plane Mach against Hall / Kliegel–Levine. Includes the Cuffel 1969 and NASA TN D-5377 experiments |
| T3 | LOX/RP-1 frozen-at-throat real gas (RocketCEA). OpenFOAM uses a JANAF cp(T) fitted to the frozen CEA expansion |
| T4 | k-ω SST + Spalding wall function: BL θ, δ* and Cf against rockmoc's integral BL and TN D-5377 surveys, BL-corrected production contour, sea-level separation (Schmucker) with a plume domain |

Not covered yet: aerospike (needs an external-flow mesher), two-phase flow, finite-rate/equilibrium chemistry in CFD.

Findings: `FINDINGS.md`. Viscous cases use a 2e-4 m first wall cell (y+ ≈ 40, Spalding wall function) to stay within an overnight budget. The cell of TN D-5377 is 2e-5 m. Separation needs the wall-resolved case `t4_sep_rao80_pc1_yp4` (2e-5 m); the coarse cases stay attached.
