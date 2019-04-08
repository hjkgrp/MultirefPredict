"""
energy_based_diag.py

Classes that calculate the energy based multireference diagnostics
"""
from abc import abstractmethod
import qcengine
import qcelemental

from MultirefPredict.io_tools import qcres_to_json, write_diagnostics_to_json
from MultirefPredict.spin import atomSpinMultDict
from MultirefPredict.cheminfo import qcelemental2OBMol
from MultirefPredict.diagnostic import Diagnostic


class EBasedDiagnostic(Diagnostic):
    def __init__(self, **kwargs):
        Diagnostic.__init__(self, **kwargs)
        self.atomized = {}

    """
    Identify unique atoms, create a dictionary for the time they repeat, 
    the qcelemental molecule instance, and the energy

    """

    def mol2atoms(self):
        print("Trying to determine the dissociation limit of the molecule\n")
        # if the dissociated atoms have been created
        if self.atomized:
            return True

        # First initialize the list of dictionary for atomized molecule
        if self.molecule.molecular_charge != 0:
            raise ValueError("Logic for assigning charge and spin states for \
                    charged molecule is still under experiment.")

        for symbol in self.molecule.symbols:
            if symbol in self.atomized:
                self.atomized[symbol]["count"] += 1
            else:
                spinDict = atomSpinMultDict()
                spinmult = spinDict.get_spinmult(symbol)
                charge = 0
                geo = [0, 0, 0]
                atom = qcelemental.models.Molecule(geometry=geo,
                                                   symbols=[symbol],
                                                   molecular_charge=charge,
                                                   molecular_multiplicity=spinmult)
                self.atomized[symbol] = {"count": 1,
                                         "molecule": atom,
                                         "energy": 0
                                         }

        # Summarize the dissociation limit
        print("Dissociation limit of the molecule:")
        print("-" * 30)
        print("{:6}".format("Atom"), "{:8}".format("Count"),
              "{:8}".format("Charge"), "{:8}".format("spin"))
        print("-" * 30)
        for symbol in self.atomized:
            print("{:^6}".format(symbol), "{:^8}".format(self.atomized[symbol]["count"]),
                  "{:^8}".format(self.atomized[symbol]["molecule"].molecular_charge),
                  "{:^8}".format(self.atomized[symbol]["molecule"].molecular_multiplicity))
        print("-" * 30)

        return True

    """
    Compute the binding energy of the molecule with given method
    """

    def energyTask(self, mol, method, program):
        task = None
        if program == "psi4":
            task = qcelemental.models.ResultInput(
                molecule=mol,
                driver="energy",
                model={"method": method, "basis": "6-31g"},
            )
        elif program == "terachem":
            tc_method = method if mol.molecular_multiplicity == 1 else "u" + method
            task = qcelemental.models.ResultInput(
                molecule=mol,
                driver="energy",
                model={"method": method, "basis": "6-31g"},
                keywords={"gpus":"1",
                          "maxit":"1500", 
                          "scf":"diis+a", 
                          "convthre":"1e-6",  
                          "precision":"double",
                          "units":"bohr",
                          "method":tc_method}
            )
        return task

    def computeBE(self, method):
        print("")
        print("Calculate atomization energy with method: ", method)

        # Caculate energy for the whole molecule
        molecule_task = self.energyTask(self.molecule, method, self.program)

        print("Evaluating the energy of the whole molecule...")
        molecule_result = qcengine.compute(molecule_task, self.program)
        if self.record:
            filename = self.rundir + "/" + self.diagnostic_type + "_" \
                       + self.molname + "_" + method + "_" + "whole" + ".json"
            qcres_to_json(molecule_result, filename=filename)
        if not molecule_result.success:
            raise RuntimeError("Quantum chemistry calculation failed.")
        molecule_energy = molecule_result.return_result
        print("Final energy of the molecule (Hartree): {:.8f}".format(molecule_energy))

        print("Evaluating the energy of the atomization limit of the molecule...")
        if not self.atomized:
            self.mol2atoms()

        # Calculate energy for each unique atom at the atomization limit
        for symbol in self.atomized:
            atom_result = None
            if self.program == "terachem" and symbol == "H":
                 atom_task = self.energyTask(self.atomized[symbol]["molecule"], method, "psi4")
                 atom_result = qcengine.compute(atom_task, "psi4")
            else:
                 atom_task = self.energyTask(self.atomized[symbol]["molecule"], method, self.program)
                 atom_result = qcengine.compute(atom_task, self.program)
            if not atom_result.success:
                raise RuntimeError("Quantum chemistry calculation failed.")
            if self.record:
                filename = self.rundir + "/" + self.diagnostic_type + "_"\
                           + self.molname + "_" + method + "_" + symbol + ".json"
                qcres_to_json(atom_result, filename=filename)
            atom_energy = atom_result.return_result
            print("Final energy of atom ", symbol, " (Hartree): {:.8f}".format(atom_energy))
            self.atomized[symbol]["energy"] = atom_energy

        # Calculate BE
        BE = 0
        for symbol in self.atomized:
            BE += self.atomized[symbol]["energy"] * self.atomized[symbol]["count"]
        print("Final energy of the atomized limit (Hartree): {:.8f}".format(BE))

        BE = BE - molecule_energy
        print("Atomization energy (Hartree): {:.8f}".format(BE))
        print("")

        return BE

    """
    Compute the diagnostic
    """

    @abstractmethod
    def computeDiagnostic(self):
        pass


class B1(EBasedDiagnostic):
    def __init__(self, **kwargs):
        EBasedDiagnostic.__init__(self, **kwargs)
        self.diagnostic_type = "B1"
        self.numBonds = qcelemental2OBMol(self.molecule).NumBonds()

    """
    Compute the B1 diagnostic
    """

    def computeDiagnostic(self):
        print("Compute B1 diagnostic of the given molecule:")
        self.molecule.pretty_print()
        self.mol2atoms()
        BE_BLYP = self.computeBE("blyp")
        BE_B1LYP = self.computeBE("b1lyp")
        print("Number of bonds in the molecule: ", self.numBonds)
        diag = (BE_BLYP - BE_B1LYP) / float(self.numBonds)
        print("\nB1 DIAGNOSTICS: {:.3f}".format(diag))
        if self.record:
            diag_dict = {self.diagnostic_type: diag}
            filename = self.rundir + "/" + self.molname + "_" + self.diagnostic_type + ".json"
            write_diagnostics_to_json(diag_dict, filename)
        return diag


class A25PBE(EBasedDiagnostic):
    def __init__(self, **kwargs):
        EBasedDiagnostic.__init__(self, **kwargs)
        self.diagnostic_type = "A25PBE"

    """
    Compute the A25PBE diagnostic
    """

    def computeDiagnostic(self):
        print("Compute A25PBE diagnostic of the given molecule:")
        self.molecule.pretty_print()
        self.mol2atoms()
        TAE_PBE = self.computeBE("pbe")
        TAE_PBE0 = self.computeBE("pbe0")
        diag = 4 * (1 - TAE_PBE0 / TAE_PBE)
        print("\nA25PBE DIAGNOSTICS: {:.3f}".format(diag))
        if self.record:
            diag_dict = {self.diagnostic_type: diag}
            filename = self.rundir + "/" + self.molname + "_" + self.diagnostic_type + ".json"
            write_diagnostics_to_json(diag_dict, filename)
        return diag


class TAE(EBasedDiagnostic):
    def __init__(self, **kwargs):
        EBasedDiagnostic.__init__(self, **kwargs)
        self.diagnostic_type = "TAE"
        if self.program != "psi4":
            raise ValueError("Support for packages other than psi4 for TAE \
                              diagnostic is to be done\n")

    """
    Compute the TAE diagnostic
    """

    def computeDiagnostic(self):
        print("Compute TAE diagnostic of the given molecule:")
        self.molecule.pretty_print()
        self.mol2atoms()
        TAE_CCSD = self.computeBE("ccsd")
        TAE_CCSDT = self.computeBE("ccsd(t)")
        diag = 100 * (1 - TAE_CCSD / TAE_CCSDT)
        print("\nTAE DIAGNOSTICS: {:.3f}".format(diag))
        if self.record:
            diag_dict = {self.diagnostic_type: diag}
            filename = self.rundir + "/" + self.molname + "_" + self.diagnostic_type + ".json"
            write_diagnostics_to_json(diag_dict, filename)
        return diag
