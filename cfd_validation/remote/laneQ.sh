#!/bin/bash
# inviscid lane
cd ~/foamstudio/nozzles/cfdval
for c in t1_rao_80_coarse t1_rao_80_fine t1_mln_M3_coarse t1_mln_M3_fine t1_cone_15_coarse t1_cone_15_fine t2_ideal_full t2_tic80_full t1_mln_g1.15 t3_perfect_contour; do rm -f $c/FAILED; bash $c/Allrun x; done
echo "$(date +%T) lane Q finished" >> queue.log
