#!/bin/bash
cd ~/foamstudio/nozzles/cfdval
until grep -q "lane A finished" queue.log; do sleep 60; done
for c in t2_ideal_full t2_tic80_full; do bash $c/Allrun x; done
echo "$(date +%T) lane C finished" >> queue.log
