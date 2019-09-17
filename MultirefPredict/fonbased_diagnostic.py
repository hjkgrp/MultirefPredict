"""
energy_based_diag.py

Classes that calculate the energy based multireference diagnostics
"""
import qcengine
import qcelemental
import math
import sys

from MultirefPredict.io_tools import qcres_to_json, write_diagnostics_to_json
from MultirefPredict.diagnostic import Diagnostic
from MultirefPredict.basis_info import molecule_to_num_AO


class FonBased(Diagnostic):
    def __init__(self, **kwargs):
        Diagnostic.__init__(self, **kwargs)
        self.diagnostic_type = "FonBased"
        self.fons = None
        self.restricted = False
        self.norb = 0
        self.ncore = 0
        self.nele = 0
        self.Hartree2K = 3.16681e-6
        if "temp_K" in kwargs.keys():
            self.temp_K = kwargs["temp_K"]
        else:
            self.temp_K = False

    def setTemperature(self, method):
        if not isinstance(method, str):
            raise ValueError("The arguement for DFT method must be of type str")
        temp = 0.01583405237
        LDAsOrGGAs = ["svwn", "pbe", "blyp", "revpbe", "bop", "b97" ]
        hybrids = {"b3lyp": 0.2, "b3lyp1": 0.2, "b3lyp5": 0.2, "b3lyp3": 0.2,\
                   "bhandhlyp": 0.5, "pbe0": 0.25, "revpbe0" : 0.25}
        # Temperature for long range correct hybrids are undefined
        if method.lower() in hybrids.keys():
            HFX = hybrids[method.lower()]
            temp = 0.06333620949 * HFX + temp
        elif method.lower() not in LDAsOrGGAs:
            raise ValueError("The temperature needed for running finite-temp DFT"\
                            +"of the given XC functional " + method + " is undefined")
        return temp

    """
    Compute the Fractional Occupation SCF
    """
    def FonTask(self, mol, program, method, basis=None):
        task = None
        #Default parameters
        basis_set = "lacvps_ecp"
        if basis is not None:
            basis_set = basis

        if not self.temp_K:
            temp = self.setTemperature(method)
            self.temp = temp
            self.temp_K = temp/self.Hartree2K
        else:
            self.temp = self.temp_K*self.Hartree2K

        print("")
        print("Fractional Occupation Number SCF with " + method + "/"+basis_set)

        print("FON temperature: {0:10.5f} (KT in atomic unit) ".format(self.temp))
        print("                 {0:10.0f} (Kelvin) ".format(self.temp_K))
        sys.stdout.flush()
        #TODO add HFX dependent temperature determination
        # Set the core orbitals as frozen and allow FON for all others
        self.norb, self.ncore = molecule_to_num_AO(mol, basis_set)
        if program == "terachem":
            tc_method = method if mol.molecular_multiplicity == 1 else "u" + method
            additional_keywords={"gpus": "1",
                "maxit": "500", 
                "scf": "diis+a", 
                "levelshift":"yes",
                "levelshiftvala":"0.25",
                "levelshiftvalb":"0.25",
                "convthre": "1e-6",  
                "precision": "double",
                "units": "bohr",
                "fon": "yes",
                "fon_method": "fermi",
                "fon_temperature": self.temp,
                "fon_print": "1",
                "method": tc_method,
                "closed": self.ncore,
                "active": self.norb-self.ncore}
            if self.wfn != None:
                if self.wfn[0]!="" and self.wfn[1]!="":
                    additional_keywords["guess"] = self.wfn[0] + " " + self.wfn[1]
                elif self.wfn[2]!="":
                    additional_keywords["guess"] = self.wfn[2]
                else:
                    print("Warning: initial guess for SCF is not provided with"\
                          + "the correct format and will be ignored")
            extra_files = {}
            if isinstance(self.extras,list):
                for extra in self.extras:
                    extra_files[extra] = ""
            task = qcelemental.models.ResultInput(
                molecule=mol,
                driver="energy",
                model={"method": tc_method, "basis": basis_set},
                keywords= additional_keywords,
                extras=extra_files
            )
        else:
            raise ValueError("FON calculation is not implemented for the requested program yet.")
        return task

    """
    Harvest FON numbers from FON calculation result
    """
    def harvestFon(self, result):
        datain = result.dict()["stdout"].split('\n')
        fons = {}
        for i in range(0,len(datain)):
            if "Total electrons:" in datain[i]:
                self.nele = int(datain[i].strip('\n').split()[2])
            if self.molecule.molecular_multiplicity == 1:
                self.restricted = True
            else:
                self.restricted = False
            if "SCF Fractional Occupations" in datain[i]:
                # Store the beta orbitals with index norb+1, norb+2, ...
                orb_offset = self.norb if len(fons) > 0 else 0
                nlines = int(math.ceil(float(self.norb)/4))
                for j in range(i+3, i+3+nlines):
                    line = datain[j].strip('\n').split()
                    for k in range(0,len(line)):
                        if k%2 == 0:
                            idx = int(line[k]) + orb_offset
                            fon = float(line[k+1])
                            if fon > 1e-6 and fon < 1.0:
                                fons[idx] = fon 
        return fons

    """
    Conduct FON calculation
    """
    def computeFon(self, method, basis=None):
        # Caculate energy for the whole molecule
        molecule_task = self.FonTask(self.molecule, self.program, method,basis)

        molecule_result = qcengine.compute(molecule_task, self.program)
        if self.record:
            filename = self.rundir + "/" + self.diagnostic_type + "_" \
                       + self.molname + "_" + method + "_" + "whole" + ".json"
            qcres_to_json(molecule_result, filename)
        if not molecule_result.success:
            raise RuntimeError("Quantum chemistry calculation failed.")
        print("FON calculation finished. Harvesting FON info.")
        self.fons = self.harvestFon(molecule_result)

    """
    compute FOD
    """
    def getFOD(self, fons, norb, nele, restricted):
        FOD = 0 
        zero_temp_homo_a = 0 
        zero_temp_homo_b = 0 
        spin_mult = self.molecule.molecular_multiplicity
        
        if restricted:
            zero_temp_homo_a = nele/2
        else:
            unpaired = spin_mult - 1 
            zero_temp_homo_b = (nele-unpaired)/2 
            zero_temp_homo_a = zero_temp_homo_b + unpaired
        
        
        factor = 2.0 if restricted else 1.0 
        
        for key,value in fons.items():
            if  key <= norb: #alpha orbitals
               f = value if key > zero_temp_homo_a else 1-value
            else: # beta orbitals
               f = value if key-norb > zero_temp_homo_b else 1-value
        
            my_FOD = f* factor
            FOD +=  my_FOD
        return FOD
    """
    compute Matito
    """
    def getMattito(self, fons, norb, restricted):
        I_ND = 0    
        I_T = 0 
        
        factor = 2.0 if restricted else 1.0 
        for key,value in fons.items():
            my_I_T = 0.25 * math.sqrt(value * (1-value)) * factor
            I_T += my_I_T
        
            my_I_ND = 0.5* value * (1-value) * factor
            I_ND +=  my_I_ND
        
        I_D = I_T - I_ND
        return I_D,I_ND
    
    def getEntanglement(self,fons, norbs, restricted):
        if not restricted:
            print("Warning: Entanglement only support calculation of natural orbital occupations.")
            print("For unrestricted entanglement calculations, Unrestricted Natural Orbitals are needed but not supported in TreaChem yet")
            return -1
        
        S=0
        for key,value in fons.items():
            S += -value*math.log(value)
        return S


    """
    Compute the diagnostic
    """

    def computeDiagnostic(self, method="PBE", basis=None):
        print("Compute FON based diagnostic of the given molecule:")
        if self.fons is None:
           self.computeFon(method, basis)
        FOD = self.getFOD(self.fons, self.norb, self.nele, self.restricted)
        I_D,I_ND = self.getMattito(self.fons, self.norb, self.restricted)
        S = self.getEntanglement(self.fons, self.norb, self.restricted)
        diag = {"FOD": FOD,
                "Mattito": {"I_D":I_D, "I_ND":I_ND},
                "Entanglement": S
               }
        if self.record:
            diag_dict = {self.diagnostic_type: diag}
            filename = self.rundir + "/" + self.molname + "_" + self.diagnostic_type + ".json"
            write_diagnostics_to_json(diag_dict, filename)
        return diag
