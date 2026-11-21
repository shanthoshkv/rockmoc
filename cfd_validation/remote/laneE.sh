#!/bin/bash
cd ~/foamstudio/nozzles/cfdval
until grep -q "lane B finished" queue.log; do sleep 60; done
for c in t3_frozen t3_equil t3_perfect_contour; do rm -f $c/FAILED $c/DONE; bash $c/Allrun x; done
echo "$(date +%T) lane E finished" >> queue.log
