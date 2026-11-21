# rockmoc CFD validation report

Generated 2026-11-21 09:11. OpenFOAM v2412 `rhoCentralFoamAxi` (Kurganov central-upwind, 2.5 deg wedge), driven by `cfd_validation/cfdval.py`.

**6 / 6 gated checks pass** (36 informational).

## Summary

| tier | case | Isp rockmoc [s] | Isp CFD [s] | dIsp | exit M (mean ± std) | MOC field M rms | Cd CFD | gated fails |
|---|---|---|---|---|---|---|---|---|
| T4 | [t4_sep_cone15_pc1](figures/t4_sep_cone15_pc1.png) | nan | 188.29 | - | 1.445 ± 1.101 | - | 1.0086 | 0 |
| T4 | [t4_sep_rao80_pc1](figures/t4_sep_rao80_pc1.png) | nan | 166.12 | - | 2.726 ± 0.807 | - | 1.0032 | 0 |
| T4 | [t4_sep_rao80_pc1.5](figures/t4_sep_rao80_pc1.5.png) | nan | 200.18 | - | 2.953 ± 0.341 | - | 1.0000 | 0 |

## T4

Viscous (k-omega SST, Spalding wall function, Sutherland viscosity): boundary-layer thicknesses, friction loss, BL-corrected contour, and sea-level flow separation with an external plume domain.

| case | quantity | rockmoc / ref | CFD | error | tol | status | note |
|---|---|---|---|---|---|---|---|
| t4_sep_cone15_pc1 | axis Mach max |err| x>0.5Rt (1-D) | 0 | 0.24161 | +0.242 | 0.02 | info |  |
| t4_sep_cone15_pc1 | axis p/p0 max |err| (1-D) | 0 | 0.73014 | +0.73 | 0.06 | info |  |
| t4_sep_cone15_pc1 | axis p/p0 rms err (1-D) | 0 | 0.37969 | +0.38 | 0.05 | info |  |
| t4_sep_cone15_pc1 | axis T/T0 max |err| (1-D) | 0 | 0.10699 | +0.107 | 0.01 | info |  |
| t4_sep_cone15_pc1 | axis T/T0 rms err (1-D) | 0 | 0.066215 | +0.0662 | 0.01 | info |  |
| t4_sep_cone15_pc1 | axis rho/rho0 max |err| (1-D) | 0 | 0.56804 | +0.568 | 0.05 | info |  |
| t4_sep_cone15_pc1 | axis rho/rho0 rms err (1-D) | 0 | 0.3071 | +0.307 | 0.04 | info |  |
| t4_sep_cone15_pc1 | axis u/a0 max |err| (1-D) | 0 | 0.20337 | +0.203 | 0.01 | info |  |
| t4_sep_cone15_pc1 | axis u/a0 rms err (1-D) | 0 | 0.063347 | +0.0633 | 0.01 | info |  |
| t4_sep_cone15_pc1 | axis Mach rms err x>0.5Rt (1-D) | 0 | 0.094187 | +0.0942 | 0.01 | info |  |
| t4_sep_cone15_pc1 | wall Mach rms err | 0 | 1 | +1 | 0.03 | info | rockmoc wall M |
| t4_sep_cone15_pc1 | separation area ratio eps_sep (Schmucker) | 4.0009 | 3.8971 | -2.59 % | 25.0 % | PASS | CFD separated |
| t4_sep_cone15_pc1 | Cf at sea level (rockmoc: Schmucker on MOC wall p, Kalt-Badal) | 1.1345 | 1.1302 | -0.38 % | 5.0 % | PASS | attached-flow theory 1.001 |
| t4_sep_cone15_pc1 | Cf at sea level vs attached-flow theory | 1.0012 | 1.1302 | +12.89 % | 5.0 % | info |  |
| t4_sep_rao80_pc1 | axis Mach max |err| x>0.5Rt (1-D) | 0 | 0.21023 | +0.21 | 0.02 | info |  |
| t4_sep_rao80_pc1 | axis p/p0 max |err| (1-D) | 0 | 0.91162 | +0.912 | 0.06 | info |  |
| t4_sep_rao80_pc1 | axis p/p0 rms err (1-D) | 0 | 0.53266 | +0.533 | 0.05 | info |  |
| t4_sep_rao80_pc1 | axis T/T0 max |err| (1-D) | 0 | 0.12337 | +0.123 | 0.01 | info |  |
| t4_sep_rao80_pc1 | axis T/T0 rms err (1-D) | 0 | 0.078243 | +0.0782 | 0.01 | info |  |
| t4_sep_rao80_pc1 | axis rho/rho0 max |err| (1-D) | 0 | 0.70169 | +0.702 | 0.05 | info |  |
| t4_sep_rao80_pc1 | axis rho/rho0 rms err (1-D) | 0 | 0.4167 | +0.417 | 0.04 | info |  |
| t4_sep_rao80_pc1 | axis u/a0 max |err| (1-D) | 0 | 0.1744 | +0.174 | 0.01 | info |  |
| t4_sep_rao80_pc1 | axis u/a0 rms err (1-D) | 0 | 0.08628 | +0.0863 | 0.01 | info |  |
| t4_sep_rao80_pc1 | axis Mach rms err x>0.5Rt (1-D) | 0 | 0.11756 | +0.118 | 0.01 | info |  |
| t4_sep_rao80_pc1 | wall Mach rms err | 0 | 1 | +1 | 0.03 | info | rockmoc wall M |
| t4_sep_rao80_pc1 | separation area ratio eps_sep (Schmucker) | 6.3395 | 5.8731 | -7.36 % | 25.0 % | PASS | CFD separated |
| t4_sep_rao80_pc1 | Cf at sea level (rockmoc: Schmucker on MOC wall p, Kalt-Badal) | 0.99546 | 0.99716 | +0.17 % | 5.0 % | PASS | attached-flow theory 0.995 |
| t4_sep_rao80_pc1 | Cf at sea level vs attached-flow theory | 0.99546 | 0.99716 | +0.17 % | 5.0 % | info |  |
| t4_sep_rao80_pc1.5 | axis Mach max |err| x>0.5Rt (1-D) | 0 | 0.21015 | +0.21 | 0.02 | info |  |
| t4_sep_rao80_pc1.5 | axis p/p0 max |err| (1-D) | 0 | 0.90996 | +0.91 | 0.06 | info |  |
| t4_sep_rao80_pc1.5 | axis p/p0 rms err (1-D) | 0 | 0.53202 | +0.532 | 0.05 | info |  |
| t4_sep_rao80_pc1.5 | axis T/T0 max |err| (1-D) | 0 | 0.12319 | +0.123 | 0.01 | info |  |
| t4_sep_rao80_pc1.5 | axis T/T0 rms err (1-D) | 0 | 0.078176 | +0.0782 | 0.01 | info |  |
| t4_sep_rao80_pc1.5 | axis rho/rho0 max |err| (1-D) | 0 | 0.70048 | +0.7 | 0.05 | info |  |
| t4_sep_rao80_pc1.5 | axis rho/rho0 rms err (1-D) | 0 | 0.41624 | +0.416 | 0.04 | info |  |
| t4_sep_rao80_pc1.5 | axis u/a0 max |err| (1-D) | 0 | 0.17463 | +0.175 | 0.01 | info |  |
| t4_sep_rao80_pc1.5 | axis u/a0 rms err (1-D) | 0 | 0.086288 | +0.0863 | 0.01 | info |  |
| t4_sep_rao80_pc1.5 | axis Mach rms err x>0.5Rt (1-D) | 0 | 0.11754 | +0.118 | 0.01 | info |  |
| t4_sep_rao80_pc1.5 | wall Mach rms err | 0 | 1 | +1 | 0.03 | info | rockmoc wall M |
| t4_sep_rao80_pc1.5 | separation area ratio eps_sep (Schmucker) | 6.3395 | 6.3395 | +0.00 % | 25.0 % | PASS | CFD: attached to lip |
| t4_sep_rao80_pc1.5 | Cf at sea level (rockmoc: Schmucker on MOC wall p, Kalt-Badal) | 1.2096 | 1.2016 | -0.66 % | 5.0 % | PASS | attached-flow theory 1.210 |
| t4_sep_rao80_pc1.5 | Cf at sea level vs attached-flow theory | 1.2096 | 1.2016 | -0.66 % | 5.0 % | info |  |

## Figures

![t4_sep_cone15_pc1](figures/t4_sep_cone15_pc1.png)
![t4_sep_rao80_pc1](figures/t4_sep_rao80_pc1.png)
![t4_sep_rao80_pc1.5](figures/t4_sep_rao80_pc1.5.png)