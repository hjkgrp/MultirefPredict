########################################################
# This script loop through structures given in xyzdir
# find the matching c0 (or ca0 cb0 for UKS) file, the 
# MO coefficient file in the directory: wfndir
# And conduct finite-temp DFT calculation
# The molden files of the FT-DFT will be saved into json output
# The results can later be pushed back into database 
#########################################################
import MultirefPredict
import qcelemental as qe
import glob
import os
import sys
import timeout_decorator

xyzdir="/home/fangliu/multiref/tmc_database/xyz/group1"
xyzlist = glob.glob(xyzdir+'/*.xyz')
timelimit = 72000
wfndir="/home/fangliu/multiref/tmc_database/wfn/group1"


#name2wfnpath = pickle.load(fwfn,fix_imports=True)

@timeout_decorator.timeout(timelimit)
def run_fod_timeout(mol, name):
    ca0_path = wfndir+"/"+name+".ca0"
    cb0_path = wfndir+"/"+name+".cb0"
    c0_path = wfndir+"/"+name+".c0"
    if os.path.exists(ca0_path) and os.path.exists(cb0_path):
        c0_path = ""
    elif os.path.exists(c0_path):
        ca0_path = ""
        cb0_path = ""
    else:
        c0_path=""
        ca0_path = ""
        cb0_path = ""
    mywfn = [ca0_path, cb0_path, c0_path]
    fon = MultirefPredict.diagnostic_factory("FOD", molecule=mol,
            program='terachem', molname=name, record=True, wfn = mywfn,
            extras=['scr/geometry.molden'])
    fon.computeDiagnostic("b3lyp")

for xyz in xyzlist:
    name = os.path.basename(xyz).strip('.xyz')
    molstr = open(xyz).read()
    mol = qe.models.Molecule.from_data(molstr, dtype='xyz+')
    print("\n\n***Calculating molecule: {}".format(name))
    sys.stdout.flush()
    try:
        run_fod_timeout(mol,name)
    except:
        print("\n!!!TeraChem timed out after {} seconds".format(timelimit)) 
    sys.stdout.flush()
