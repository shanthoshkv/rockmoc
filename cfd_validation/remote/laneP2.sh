#!/bin/bash
cd ~/foamstudio/nozzles/cfdval
for c in t4_prod_bl t4_mln_inviscid_wall_visc t4_sep_rao80_pc1 t4_sep_cone15_pc1 t4_sep_rao80_pc1.5 t4_d5377_adiabatic; do rm -f $c/FAILED; bash $c/Allrun x; done
echo "$(date +%T) lane P2 finished" >> queue.log
