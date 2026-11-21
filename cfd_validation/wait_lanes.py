import sys, time
sys.path.insert(0, r"C:\Users\kvsha\Documents\Learning\openfoam scripting website\scripts")
import remote_run as rr
lanes = sys.argv[1:]
while True:
    try:
        q = rr.wsl("cat ~/foamstudio/nozzles/cfdval/queue.log")
    except SystemExit:
        q = ""
    if all(f"lane {l} finished" in q for l in lanes):
        print("finished", lanes); break
    time.sleep(120)
