import pytest
import qcelemental
import qcengine
from qcelemental.testing import compare_recursive
from MultirefPredict.casbased_diagnostic import obtainC0fromstr


def test_casscf():
    mol = qcelemental.models.Molecule(geometry=[0, 0, 0],
                                      symbols=["O"],
                                      molecular_charge=0,
                                      molecular_multiplicity=3)
    casscf_task = qcelemental.models.ResultInput(molecule=mol,
                                                 driver="energy",
                                                 model={"method": "casscf", "basis": "cc-pvdz"},
                                                 keywords={"reference": "rohf"}
                                                 )
    results = qcengine.compute(casscf_task, "psi4")
    energy = results.return_result
    expected = -74.91174384584482
    Thre = 1e-6
    assert compare_recursive(energy, expected, atol=Thre)


def test_C0():
    mol = qcelemental.models.Molecule(geometry=[[0, 0, 0], [0.7, 0, 0]],
                                      symbols=["H", "H"],
                                      molecular_charge=0,
                                      molecular_multiplicity=1)
    casscf_task = qcelemental.models.ResultInput(molecule=mol,
                                                 driver="energy",
                                                 model={"method": "casscf", "basis": "cc-pvdz"},
                                                 keywords={"reference": "rohf"}
                                                 )
    results = qcengine.compute(casscf_task, "psi4")
    outfile = results.stdout.split("\n")
    diag = obtainC0fromstr(outfile)
    expected = 0.995293
    Thre = 1e-3
    assert compare_recursive(abs(diag), expected, atol=Thre)
