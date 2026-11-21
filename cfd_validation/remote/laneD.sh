#!/bin/bash
cd ~/foamstudio/nozzles/cfdval
until grep -q "lane C finished" queue.log; do sleep 60; done
for c in t1_mln_g1.15; do bash $c/Allrun x; done
echo "$(date +%T) lane D finished" >> queue.log
