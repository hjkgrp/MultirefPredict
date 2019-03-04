"""
MultirefPredict.py
Automated workflow to predict multireference character of molecules

Handles the primary functions
"""
from abc import ABC, abstractmethod


def diagostic_factory(xyz, diagnostic_type, program):
    """
    factory class for calculating diagnostics

    Parameters
    ----------
    xyz : str
        Path to the input xyz file for the molecule to be studied
    diagnostic_type: {"B1","A25PBE","TAE"}
        The type of diagnostic to calculate
    program: {"psi4"} 
        The program to calculate the diagnostic

    Returns
    -------
    cls_instance: instance of a class for calculate specific diagnostics
    """
    cls_dict = dict(B1=B1,A25TPBE=A25PBE,TAE=TAE)
    
    if diagnostic_type not in cls_dict.keys():
        raise Exception("Diagnostic type not found")

    cls = cls_dict[diagnostic_type]
    cls_instance = cls(xyz,program)
    return cls_instance

class Diagnostic(ABC):
    @abstractmethod
    def compute_diagnostic(self):
        pass

class B1(Diagnostic):
    def __init__(self, xyz, program):
        self.program = program
        self.xyz = xyz
        #TODO: fill in

    def compute_diagnostic(self):
        pass

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print("Welcome to Multireference Diagnostics")
