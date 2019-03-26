"""
multirefpredict.py
Automated workflow to predict multireference character of molecules in quantum chemistry calculation

Handles the primary functions
"""
from abc import ABC, abstractmethod
import qcelemental
from MultirefPredict.cheminfo import xyzfile2qcelemental
from MultirefPredict.io_tools import ensure_dir
from MultirefPredict.globs import keywords_avail, available_programs


class Diagnostic(ABC):
    ### Attributes that are in common for different diagnostics should go in here.
    def __init__(self, **kwargs):
        print("inputs dictionary:", kwargs)
        self.keywords_avail = keywords_avail
        self.diagnostic_type = kwargs["diagnostic_type"]
        self.xyzfile = False
        self.molname = kwargs["molname"] if "molname" in kwargs.keys() else "undef"
        self.rundir = kwargs["rundir"] if "rundir" in kwargs.keys() else "./"
        self.record = kwargs["record"] if "record" in kwargs.keys() else False
        self.program = kwargs["program"] if "program" in kwargs.keys() else "psi4"
        self.initialization_check(**kwargs)

    def initialization_check(self, **kwargs):
        ensure_dir(self.rundir)
        for key, value in kwargs.items():
            if key not in self.keywords_avail:
                raise KeyError("Energy based diagnostic: unrecoganized key")
        if self.program not in available_programs:
            raise ValueError("Energy based diagnostic: specified program is not supported yet")
        if "xyzfile" in kwargs.keys():
            self.xyzfile = kwargs["xyzfile"]
            self.charge = kwargs["charge"] if "charge" in kwargs.keys() else 0
            self.spinmult = kwargs["spinmult"] if "spinmult" in kwargs.keys() else 1
            self.molecule = xyzfile2qcelemental(self.xyzfile, charge=self.charge, spinmult=self.spinmult)
        elif "molecule" in kwargs.keys():
            self.molecule = kwargs['molecule']
            if not isinstance(self.molecule, qcelemental.models.Molecule):
                raise TypeError("Argument molecule must be a molecule instance")
            if self.xyzfile:
                raise ValueError(
                    "Cannot have both the xyzfile and 'qcelemental.models.Molecule' as inputs the same time.")
        else:
            raise TypeError("Must input either a xyzfile file or an instance of 'qcelemental.models.Molecule'")

    @abstractmethod
    def computeDiagnostic(self):
        pass
