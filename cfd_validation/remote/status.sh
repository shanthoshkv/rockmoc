#!/bin/bash
cd ~/foamstudio/nozzles/cfdval
uptime
tail -30 queue.log
echo "done: $(ls -d */DONE 2>/dev/null | wc -l)"
ls -d */FAILED 2>/dev/null
for f in */log.solver; do
  d=$(dirname "$f")
  [ -f "$d/DONE" ] && continue
  t=$(grep '^Time =' "$f" | tail -1)
  e=$(grep -o 'endTime [^;]*' "$d/system/controlDict")
  echo "$d $t / $e"
done
