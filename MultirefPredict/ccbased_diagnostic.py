"""
energy_based_diag.py

Classes that calculate the energy based multireference diagnostics
"""
import qcengine
import qcelemental
from MultirefPredict.diagnostic import Diagnostic
from MultirefPredict.io_tools import qcres_to_json, write_diagnostics_to_json


class CCBased(Diagnostic):
    def __init__(self, **kwargs):
        Diagnostic.__init__(self, **kwargs)
        self.diagnostic_type = "CCbased"
        self.result = False
        if self.program != "psi4":
            raise ValueError("Support for packages other than psi4 for CCBased \
                              diagnostics is to be done\n")

    """
    Do CCSD(T) calculation for the given molecule
    Input: self
    Ouput: QCElemental Result
    """

    def computeCCSDT(self):
        print("")
        print("Preparing CCSD(T) calculation")
        method = "ccsd(t)"
        basis = "cc-pvdz"
        if self.program != "psi4":
            raise ValueError("Support for packages other than psi4 is to be done\n")

        # Caculate energy for the whole molecule
        molecule_task = qcelemental.models.ResultInput(
            molecule=self.molecule,
            driver="energy",
            model={"method": method, "basis": basis},
        )
        print("Evaluating the energy of the whole molecule...")
        self.result = qcengine.compute(molecule_task, "psi4")
        print("self.result:", self.result)
        if not self.result.success:
            raise RuntimeError("Quantum chemistry calculation failed.")
        if self.record:
            filename = self.rundir + "/" + self.diagnostic_type + "_" + self.molname + "_" + "ccsd(t)" + "_" + "whole" + ".json"
            qcres_to_json(self.result, filename=filename)

    """
    Compute the CCBased diagnostic
    """

    def computeDiagnostic(self):
        print("Compute CC based diagnostics of the given molecule:")
        self.molecule.pretty_print()
        self.computeCCSDT()

        if not self.result.success:
            raise RuntimeError("Quantum chemistry calculation failed.")

        T1 = self.result.extras['local_qcvars']['CC T1 DIAGNOSTIC']
        D1 = self.result.extras['local_qcvars']['CC D1 DIAGNOSTIC']
        D2 = self.result.extras['local_qcvars']['CC D2 DIAGNOSTIC']
        NewD1 = self.result.extras['local_qcvars']['CC NEW D1 DIAGNOSTIC']
        diag = {"T1": T1, "D1": D1, "D2": D2, "New D1": NewD1}
        print("\nCCBased DIAGNOSTICS:", diag)
        if self.record:
            filename = self.rundir + "/" + self.molname + "_" + self.diagnostic_type + ".json"
            write_diagnostics_to_json(diag, filename)
        return diag
