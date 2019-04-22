"""
Unit and regression test for the MultirefPredict package.
"""

# Import package, test suite, and other packages as needed
import MultirefPredict
import pytest
import sys
import os
from .compare import fuzzyEqual

def test_MultirefPredict_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "MultirefPredict" in sys.modules

def test_diagnostic_factory(qcelemental_water):
    calculator = MultirefPredict.diagnostic_factory("B1", molecule=qcelemental_water)
    diag = calculator.computeDiagnostic()
    B1Thre = 1e-6
    expected = 0.006334860365228678
    assert fuzzyEqual(diag, expected, B1Thre)

    calculator = MultirefPredict.diagnostic_factory("A25PBE", molecule=qcelemental_water)
    diag = calculator.computeDiagnostic()
    A25PBEThre = 1e-6
    expected = 0.1626572016077259
    assert fuzzyEqual(diag, expected, A25PBEThre)

    calculator = MultirefPredict.diagnostic_factory("TAE", molecule=qcelemental_water)
    diag = calculator.computeDiagnostic()
    TAEThre = 1e-6
    expected = 0.28078682517214126
    assert fuzzyEqual(diag, expected, TAEThre)

    #Thre = 1e-6
    #calculator = MultirefPredict.diagnostic_factory("CCBased", molecule=qcelemental_water)
    #diag = calculator.computeDiagnostic()
    #expected = {'T1': 0.005903453140867589, 'D1': 0.012965638464472996, 'D2': 0.12885433260716933, 'New D1': 0.012965638464472996}
    #print(diag)
    #for key,value in diag.items():
    #    assert fuzzyEqual(value, expected[key], Thre)

    #calculator = MultirefPredict.diagnostic_factory("T1", molecule=qcelemental_water)
    #diag = calculator.computeDiagnostic()
    #expected = {'T1': 0.005903453140867589, 'D1': 0.012965638464472996, 'D2': 0.12885433260716933, 'New D1': 0.012965638464472996}
    #print(diag)
    #for key,value in diag.items():
    #    assert fuzzyEqual(value, expected[key], Thre)

    #calculator = MultirefPredict.diagnostic_factory("D1", molecule=qcelemental_water)
    #diag = calculator.computeDiagnostic()
    #expected = {'T1': 0.005903453140867589, 'D1': 0.012965638464472996, 'D2': 0.12885433260716933, 'New D1': 0.012965638464472996}
    #print(diag)
    #for key,value in diag.items():
    #    assert fuzzyEqual(value, expected[key], Thre)

    #calculator = MultirefPredict.diagnostic_factory("D2", molecule=qcelemental_water)
    #diag = calculator.computeDiagnostic()
    #expected = {'T1': 0.005903453140867589, 'D1': 0.012965638464472996, 'D2': 0.12885433260716933, 'New D1': 0.012965638464472996}
    #print(diag)
    #for key,value in diag.items():
    #    assert fuzzyEqual(value, expected[key], Thre)

