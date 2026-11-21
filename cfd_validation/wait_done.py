import sys, time
sys.path.insert(0, r"C:\Users\kvsha\Documents\Learning\openfoam scripting website\scripts")
import remote_run as rr
target = int(sys.argv[1])
while True:
    try:
        n = int(rr.wsl("ls -d ~/foamstudio/nozzles/cfdval/*/DONE 2>/dev/null | wc -l").strip() or 0)
    except SystemExit:
        n = -1
    if n >= target:
        print("done count", n); break
    time.sleep(60)
