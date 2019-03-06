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

available_programs =  ["psi4"]
class EBasedDiagnostic(Diagnostic):
    def __init__(self, **kwargs):

        for key,value in kwargs.items():
            if key not in ["molecule","program"]:
                raise KeyError("Energy based diagnostic: unrecoganized key")

        self.molecule = kwargs['molecule']

        if not isinstance(self.molecule,qcelemental.models.Molecule):
            raise TypeError("Argument molecule must be a molecule instance")

        if "program" in kwargs:
            self.program = kwargs["program"]
        else:
            self.program = "psi4"

        if self.program not in available_programs:
            raise ValueError("Energy based diagnostic: specified program is not supported yet")

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
                        molecular_charge = charge, 
                        molecular_multiplicity = spinmult )
                self.atomized[symbol] = { "count": 1,
                                          "molecule": atom,
                                          "energy": 0
                                        }

        # Summarize the dissociation limit
        print("Dissociation limit of the molecule:")
        print("-"*30)
        print("{:6}".format("Atom"),"{:8}".format("Count"),
              "{:8}".format("Charge"),"{:8}".format("spin"))
        print("-"*30)
        for symbol in self.atomized: 
            print("{:^6}".format(symbol),"{:^8}".format(self.atomized[symbol]["count"]),
                  "{:^8}".format(self.atomized[symbol]["molecule"].molecular_charge),
                  "{:^8}".format(self.atomized[symbol]["molecule"].molecular_multiplicity))
        print("-"*30)

        return True

    """
    Compute the binding energy of the molecule with given method
    """
    def computeBE(self, method):
        print("")
        print("Calculate atomization energy with method: ", method)
        if self.program != "psi4":
            raise ValueError("Support for packages other than psi4 is to be done\n")

        # Caculate energy for the whole molecule
        molecule_task = qcelemental.models.ResultInput (
                molecule = self.molecule,
                driver = "energy",
                model = {"method" : method, "basis" : "6-31g" },
        )

        print("Evaluating the energy of the whole molecule...")
        molecule_result = qcengine.compute(molecule_task, "psi4")
        if not molecule_result.success:
            raise RuntimeError("Quantum chemistry calculation failed.")
        molecule_energy = molecule_result.return_result
        print("Final energy of the molecule (Hartree): {:.8f}".format(molecule_energy))

        print("Evaluating the energy of the atomization limit of the molecule...")
        if not self.atomized:
            self.mol2atoms()

        # Calculate energy for each unique atom at the atomization limit
        for symbol in self.atomized:
            atom_task = qcelemental.models.ResultInput (
                molecule = self.atomized[symbol]["molecule"],
                driver = "energy",
                model = {"method" :  method, "basis" :  "6-31g"},
            )
            atom_result = qcengine.compute(atom_task, "psi4")
            if not atom_result.success:
                raise RuntimeError("Quantum chemistry calculation failed.")
            atom_energy = atom_result.return_result
            print("Final energy of atom ",symbol," (Hartree): {:.8f}".format(atom_energy))
            self.atomized[symbol]["energy"] = atom_energy

        #Calculate BE
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
        return diag

class A25PBE(EBasedDiagnostic):
    def __init__(self, **kwargs):
        EBasedDiagnostic.__init__(self, **kwargs)

    """
    Compute the A25PBE diagnostic
    """
    def computeDiagnostic(self):
        print("Compute A25PBE diagnostic of the given molecule:")
        self.molecule.pretty_print()
        self.mol2atoms()
        TAE_PBE = self.computeBE("pbe")
        TAE_PBE0 = self.computeBE("pbe0")
        diag = 4*(1-TAE_PBE0/TAE_PBE)
        print("\nA25PBE DIAGNOSTICS: {:.3f}".format(diag))
        return diag

class TAE(EBasedDiagnostic):
    def __init__(self, **kwargs):
        EBasedDiagnostic.__init__(self, **kwargs)

    """
    Compute the TAE diagnostic
    """
    def computeDiagnostic(self):
        print("Compute TAE diagnostic of the given molecule:")
        self.molecule.pretty_print()
        self.mol2atoms()
        TAE_CCSD = self.computeBE("ccsd")
        TAE_CCSDT = self.computeBE("ccsd(t)")
        diag = 100*(1-TAE_CCSD/TAE_CCSDT)
        print("\nTAE DIAGNOSTICS: {:.3f}".format(diag))
        return diag
