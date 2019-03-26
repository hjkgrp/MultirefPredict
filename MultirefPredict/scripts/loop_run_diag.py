import os
import shutil
from MultirefPredict.io_tools import submit_job, get_jobids, ensure_dir

numjobs = len(get_jobids(req=False))
maxjobs = 150
print("Number of live jobs in queue:", numjobs)
cwd = os.getcwd()
for dirpath, dirs, files in os.walk(cwd):
    if dirpath == cwd:
        for file in files:
            if file.split(".")[-1] == "xyz":
                basepath = cwd + '/DFT_based/' + ".".join(file.split(".")[:-1])
                ensure_dir(basepath)
                shutil.copy("jobscript_sge.sh", basepath + "/jobscript.sh")
                shutil.copy("run_diagnostics_ebased.py", basepath + "/run_diagnostics_ebased.py")
                shutil.copy(dirpath + '/' + file, basepath + '/geo.xyz')
                submit_job(basepath, cwd)
                numjobs += 1
            if numjobs > maxjobs:
                print("Maximum number of jobs is reached.")
                quit()
print("All jobs are submitted.")
