#!/bin/bash
cd ~/foamstudio/nozzles/cfdval
for c in t1_rao_60 t1_rao_80 t1_rao_100 t1_raotoc_80 t1_cone_15 t1_cone_10 t1_ideal_M3 t1_tic_60 t1_tic_80 t1_mln_M3 t1_mln_M2 t1_mln_M4 t1_mln_M5 t1_mln_g1.15 t1_mln_g1.4 t1_mln_n10 t1_mln_n20 t1_mln_n80 t1_prod_mln t1_raoctl_80 t1_planar_M2.4 t2_rao80_rc1 t2_rao80_rc0.625 t2_rao80_rc2 t2_ideal_full t2_tic80_full t2_mln_full t2_cuffel t2_d5377_inv t3_frozen t3_equil t3_perfect_contour t1_mln_M3_coarse t1_mln_M3_fine t1_rao_80_coarse t1_rao_80_fine t1_cone_15_coarse t1_cone_15_fine t4_d5377_cooled; do rm -f $c/FAILED; bash $c/Allrun x; done
echo "$(date +%T) lane Q2 finished" >> queue.log
