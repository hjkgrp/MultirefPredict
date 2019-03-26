"""
casscf_based_diag.py

Classes that calculate the casscf based multireference diagnostics
"""
from abc import abstractmethod
import numpy as np
import qcengine
import qcelemental
from MultirefPredict.diagnostic import Diagnostic
from MultirefPredict.io_tools import qcres_to_json, write_diagnostics_to_json

available_programs = ["psi4"]


def obtainC0fromstr(outfile):
    C0 = 'undef'
    for line in outfile:
        if "*   1" in line:
            try:
                C0 = float(line.split()[2])
            except:
                C0 = np.nan
            break
    return C0


class CASBasedDiagnostic(Diagnostic):
    def __init__(self, **kwargs):
        Diagnostic.__init__(self, **kwargs)
        self.results = False
        self.C0 = np.nan

    def compute(self):
        print(self.molecule)
        casscf_task = qcelemental.models.ResultInput(molecule=self.molecule,
                                                     driver="energy",
                                                     model={"method": "casscf", "basis": "cc-pvdz"},
                                                     keywords={"reference": "rohf"}
                                                     )
        self.results = qcengine.compute(casscf_task, "psi4")
        if not self.results.success:
            raise RuntimeError("Quantum chemistry calculation failed.")
        else:
            self.get_C0()
        if self.record:
            filename = self.rundir + "/" + self.diagnostic_type + "_" + self.molname + "_" + "casscf" + "_" + "whole" + ".json"
            qcres_to_json(self.results, filename=filename)

    def get_C0(self):
        outfile = self.results.stdout.split("\n")
        self.C0 = obtainC0fromstr(outfile)

    """
    Compute the diagnostic
    """

    @abstractmethod
    def computeDiagnostic(self):
        pass


class C0(CASBasedDiagnostic):
    def __init__(self, **kwargs):
        CASBasedDiagnostic.__init__(self, **kwargs)
        self.diagnostic_type = "C0"

    """
    Compute the B1 diagnostic
    """

    def computeDiagnostic(self):
        print("Compute C0 diagnostic of the given molecule:")
        self.molecule.pretty_print()
        self.compute()
        diag = dict()
        diag.update({"C0": round(self.C0, 6),
                     "C0^2": round(self.C0 * self.C0, 6)
                     })
        print("\nC0 DIAGNOSTICS: ", diag)
        if self.record:
            filename = self.rundir + "/" + self.molname + "_" + self.diagnostic_type + ".json"
            write_diagnostics_to_json(diag, filename)
        return diag
