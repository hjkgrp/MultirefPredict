"""
Unit and regression test for the MultirefPredict package.
"""

# Import package, test suite, and other packages as needed
import MultirefPredict
import pytest
import sys
import os
from .compare import fuzzyEqual

@pytest.mark.skipif("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true", reason="Skipping this test on Travis CI.")
def test_MultirefPredict_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "MultirefPredict" in sys.modules

@pytest.mark.skipif("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true", reason="Skipping this test on Travis CI.")
def test_diagnostic_factory(qcelemental_water):
    calculator = MultirefPredict.diagnostic_factory("B1", molecule=qcelemental_water)
    result = calculator.computeDiagnostic()
    B1Thre = 1e-6
    assert fuzzyEqual(result, 0.006334860365228678, B1Thre)
