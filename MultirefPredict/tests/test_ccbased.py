"""
Unit and regression test for the MultirefPredict package.
"""

# Import package, test suite, and other packages as needed
import pytest
from MultirefPredict.ccbased_diagnostic import CCBased
from .compare import fuzzyEqual


@pytest.fixture(scope="class")
def cc_water(qcelemental_water):
    ccbased = CCBased(molecule=qcelemental_water)
    return ccbased

@pytest.fixture(scope="class")
def cc_water_not_psi4(qcelemental_water):
    ccbased = CCBased(molecule=qcelemental_water)
    ccbased.program = "terachem"
    return ccbased

def test_cc_init(qcelemental_water, xyz_water):
    ccbased_calculator = CCBased(molecule=qcelemental_water)
    assert ccbased_calculator.molecule == qcelemental_water
    assert ccbased_calculator.program == "psi4"

    with pytest.raises(TypeError):
         ccbased_calculator = CCBased(molecule="NothingHere")

    with pytest.raises(KeyError):
         ccbased_calculator = CCBased()

    with pytest.raises(ValueError):
         ccbased_calculator = CCBased(molecule = qcelemental_water, 
                                      program = "terachem")
def test_cc_computeCCSDT(cc_water,cc_water_not_psi4):
    cc_water.computeCCSDT()
    assert cc_water.result.success
    
    with pytest.raises(ValueError):
        cc_water_not_psi4.computeCCSDT()

def test_cc_computeDiagnostic(cc_water):
    Thre = 1e-6
    diag = cc_water.computeDiagnostic()
    expected = {'T1': 0.005903453140867589, 'D1': 0.012965638464472996, 'D2': 0.12885433260716933, 'New D1': 0.012965638464472996}
    print(diag)
    for key,value in diag.items():
        assert fuzzyEqual(value, expected[key], Thre)

