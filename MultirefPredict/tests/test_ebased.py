"""
Unit and regression test for the MultirefPredict package.
"""

# Import package, test suite, and other packages as needed
import pytest
import qcelemental
from MultirefPredict.ebased_diagnostic import B1, A25PBE, TAE
from qcengine.testing import using_terachem
from qcelemental.testing import compare_recursive
import sys
import qcengine
import psi4

def test_b1_init(qcelemental_water, xyz_water):
    b1_calculator = B1(molecule = qcelemental_water)
    assert b1_calculator.molecule == qcelemental_water
    assert b1_calculator.numBonds == 2
    assert b1_calculator.program == "psi4"
    assert b1_calculator.atomized == {}

    with pytest.raises(TypeError):
         b1_calculator = B1(molecule="NothingHere")

    with pytest.raises(KeyError):
         b1_calculator = B1()


@pytest.fixture(scope="class")
def b1_water(qcelemental_water):
    b1 = B1(molecule=qcelemental_water)
    return b1


def test_b1_mol2atoms(b1_water):
    b1_water.mol2atoms()
    assert len(b1_water.atomized) == 2
    assert "O" in b1_water.atomized
    assert b1_water.atomized["O"]["count"] == 1
    assert b1_water.atomized["O"]["molecule"].symbols == ["O"]
    assert b1_water.atomized["O"]["molecule"].molecular_charge == 0
    assert b1_water.atomized["O"]["molecule"].molecular_multiplicity == 3

    assert "H" in b1_water.atomized
    assert b1_water.atomized["H"]["count"] == 2
    assert b1_water.atomized["H"]["molecule"].symbols == ["H"]
    assert b1_water.atomized["H"]["molecule"].molecular_charge == 0
    assert b1_water.atomized["H"]["molecule"].molecular_multiplicity == 2

def test_qcengine():
    assert "qcengine" in sys.modules

def test_psi4():
    assert "psi4" in sys.modules
    psi4.set_memory('500 MB')

    h2o = psi4.geometry("""
            O
            H 1 0.96
            H 1 0.96 2 104.5
            """)

    res=psi4.energy('scf/cc-pvdz')
    assert res!=0

def test_energy_prep(b1_water):
    # Caculate energy for the whole molecule
    dft_functional = "blyp"
    molecule_task = {
            "schema_name": "qcschema_input",
            "schema_version": 1,
            "molecule": b1_water.molecule,
            "driver": "energy",
            "model": {"method": dft_functional, "basis": "6-31g" 
                     },
    }

    print("Evaluating the energy of the whole molecule...")
    molecule_result = qcengine.compute(molecule_task, "psi4")
    assert True

def test_b1_computeBE(b1_water):
    EnergyThre = 1e-6

    BE_blyp = b1_water.computeBE("blyp")
    expected = 0.33100268529368293
    assert compare_recursive(BE_blyp, expected, atol = EnergyThre)

    BE_b1lyp = b1_water.computeBE("b1lyp")
    expected = 0.3183329645632256
    assert compare_recursive(BE_b1lyp, expected, atol = EnergyThre)

def test_b1_computeDiagnostic(b1_water):
    B1Thre = 1e-6
    diag = b1_water.computeDiagnostic()
    expected = 0.006334860365228678
    assert compare_recursive(diag, expected, atol = B1Thre)

#@pytest.fixture(scope="class")
#def b1_water_terachem(qcelemental_water):
#    b1 = B1(molecule=qcelemental_water, program="terachem")
#    return b1
#
#def test_b1_terachem(b1_water_terachem):
#    B1Thre = 1e-6
#    diag = b1_water_terachem.computeDiagnostic()
#    expected = 0.006334860365228678
#    assert compare_recursive(diag, expected, atol = B1Thre)

@pytest.fixture(scope="class")
def a25pbe_water(qcelemental_water):
    a25pbe = A25PBE(molecule=qcelemental_water)
    return a25pbe

def test_a25pbe_computeDiagnostic(a25pbe_water):
    A25PBEThre = 1e-3
    diag = a25pbe_water.computeDiagnostic()
    expected = 0.1626572016077259
    assert compare_recursive(diag, expected, atol = A25PBEThre)

@pytest.fixture(scope="class")
def a25_water_terachem(qcelemental_water):
    a25 = A25PBE(molecule=qcelemental_water, program="terachem")
    return a25

@using_terachem
def test_a25_terachem(a25_water_terachem):
    A25PBEThre = 1e-3
    diag = a25_water_terachem.computeDiagnostic()
    expected = 0.1626572016077259
    assert compare_recursive(diag, expected, atol = A25PBEThre)

@pytest.fixture(scope="class")
def tae_water(qcelemental_water):
    tae = TAE(molecule=qcelemental_water)
    return tae 

def test_a25tae_computeDiagnostic(tae_water):
    TAEThre = 1e-3
    diag = tae_water.computeDiagnostic()
    expected = 0.281577634
    assert compare_recursive(diag, expected, atol = TAEThre)

