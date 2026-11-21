cd ~/foamstudio/nozzles/cfdval
for d in */; do [ -f $d/Allrun ] && cp /mnt/c/Users/SHANTHOSH/xfer/Allrun $d/Allrun; done
echo fixed
