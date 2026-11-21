#!/bin/bash
cd ~/foamstudio/nozzles/cfdval
for c in t4_sep_rao80_pc1_yp4; do rm -f $c/FAILED; bash $c/Allrun x; done
echo "$(date +%T) lane R finished" >> queue.log
