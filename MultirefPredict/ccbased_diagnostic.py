"""
energy_based_diag.py

Classes that calculate the energy based multireference diagnostics
"""
from abc import abstractmethod
import qcengine
import qcelemental
from MultirefPredict.spin import atomSpinMultDict
from MultirefPredict.cheminfo import qcelemental2OBMol
from MultirefPredict.diagnostic import Diagnostic


class CCBased(Diagnostic):
    def __init__(self, **kwargs):
        Diagnostic.__init__(self, **kwargs)
        self.diagnostic_type = "CCbased"
        self.result = False

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
        molecule_task = qcelemental.models.ResultInput (
                molecule = self.molecule,
                driver = "energy",
                model = {"method" : method, "basis" : basis },
        )

        print("Evaluating the energy of the whole molecule...")
        molecule_result = qcengine.compute(molecule_task, "psi4")
        return molecule_result

    """
    Compute the B1 diagnostic
    """
    def computeDiagnostic(self):
        print("Compute CC based diagnostics of the given molecule:")
        self.molecule.pretty_print()
        self.result = self.computeCCSDT()

        if not self.result.success:
            raise RuntimeError("Quantum chemistry calculation failed.")

        T1=self.result.extras['local_qcvars']['CC T1 DIAGNOSTIC']
        D1=self.result.extras['local_qcvars']['CC D1 DIAGNOSTIC']
        D2=self.result.extras['local_qcvars']['CC D2 DIAGNOSTIC']
        NewD1=self.result.extras['local_qcvars']['CC NEW D1 DIAGNOSTIC']
        diag = {"T1":T1, "D1":D1, "D2":D2, "New D1":D1}
        return diag
