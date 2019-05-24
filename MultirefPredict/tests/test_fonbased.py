"""
Unit and regression test for the MultirefPredict package.
"""

# Import package, test suite, and other packages as needed
import pytest
from MultirefPredict.fonbased_diagnostic import FonBased
from qcengine.testing import using_terachem
from qcelemental.testing import compare_recursive


@pytest.fixture(scope="class")
def fon_cu_complex(qcelemental_cu_complex):
    fonbased = FonBased(molecule=qcelemental_cu_complex, program="terachem")
    return fonbased

@using_terachem
def test_fon_computeFon(fon_cu_complex):
    result = fon_cu_complex.computeFon("PBE")
    assert fon_cu_complex.fons is not None

@using_terachem
def test_fon_computeDiagnostic(fon_cu_complex):
    Thre = 1e-3
    diag = fon_cu_complex.computeDiagnostic()
    expected = {'FOD': 0.980722,
                'Mattito':{ 'I_D': 0.473731315028, 'I_ND':  0.431365088635}, 
                'Entanglement': 0.732773215873}
    print(diag)
    assert compare_recursive(diag, expected, atol=Thre)

