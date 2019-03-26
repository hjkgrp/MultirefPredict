import os
import subprocess
import time
import json


def qcres_to_json(res, filename):
    outdict = res.dict()
    outdict['molecule']['geometry'] = outdict['molecule']['geometry'].flatten().tolist()
    f = open(filename, 'w')
    json.dump(outdict, f)
    f.close()


def write_diagnostics_to_json(diag_dict, filename):
    f = open(filename, 'w')
    json.dump(diag_dict, f)
    f.close()


def ensure_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


def submit_job(rundir, cwd):
    print("Submitting jobs: %s" % rundir)
    os.chdir(rundir)
    run_cmd = 'qsub jobscript.sh'
    q = subprocess.Popen(run_cmd, shell=True, stdout=subprocess.PIPE)
    ll = q.communicate()[0].decode("utf-8")
    print('ll:', ll)
    os.chdir(cwd)
    time.sleep(1)


def get_jobids(req):
    run_cmd = 'qstat'
    q = subprocess.Popen(run_cmd, shell=True, stdout=subprocess.PIPE)
    txt = q.communicate()[0].decode("utf-8")
    jobs = txt.split("\n")
    joblist = []
    for job in jobs[2:-1]:
        # if str(job.split()[4]) == "r":
        joblist.append([str(job.split()[0]), str(job.split()[7])])
    return joblist
