"""
Unit and regression test for the MultirefPredict package.
"""

# Import package, test suite, and other packages as needed
import pytest
import qcelemental
from MultirefPredict.b1 import B1
from .compare import fuzzyEqual

#def test_b1():
#    assert "B1" in dir(MultirefPredict.b1)

def test_b1_init(qcelemental_water, xyz_water):
    b1_calculator = B1(molecule = qcelemental_water)
    assert b1_calculator.molecule == qcelemental_water
    assert b1_calculator.numBonds == 2
    assert b1_calculator.program == "psi4"
    assert b1_calculator.atomized == {}

    with pytest.raises(TypeError):
         b1_calculator = B1(molecule=xyz_water)

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

def test_b1_computeBE(b1_water):
    EnergyThre = 1e-6

    BE_blyp = b1_water.computeBE("blyp")
    assert fuzzyEqual(BE_blyp,0.33100268529368293, EnergyThre)

    BE_b1lyp = b1_water.computeBE("b1lyp")
    assert fuzzyEqual(BE_b1lyp,0.3183329645632256,EnergyThre)

def test_b1_computeDiagnostic(b1_water):
    B1Thre = 1e-6
    B1 = b1_water.computeDiagnostic()
    assert fuzzyEqual(B1,0.006334860365228678, B1Thre)
