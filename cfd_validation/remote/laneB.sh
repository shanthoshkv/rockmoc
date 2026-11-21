#!/bin/bash
cd ~/foamstudio/nozzles/cfdval
for c in t1_cone_15 t1_ideal_M3 t1_mln_g1.4 t1_mln_M4 t1_mln_n10 t1_mln_n20 t1_planar_M2.4 t1_rao_100 t1_rao_60 t1_raoctl_80 t1_tic_60 t1_tic_80 t2_d5377_inv t2_rao80_rc1 t3_equil t3_perfect_contour t4_d5377_cooled t4_prod_bl t4_sep_rao80_pc1 t4_sep_rao80_pc1.5 t1_rao_80_coarse t1_rao_80_fine; do bash $c/Allrun x; done
echo "$(date +%T) lane B finished" >> queue.log
