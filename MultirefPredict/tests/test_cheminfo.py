"""
Unit and regression test for the MultirefPredict package.
"""

# Import package, test suite, and other packages as needed
import pytest
import sys
import qcelemental
from MultirefPredict import cheminfo


def test_qcelemental2Xyz(qcelemental_water, xyz_water):
    assert cheminfo.qcelemental2Xyz(qcelemental_water) == xyz_water

def test_qcelemental2OBMol(qcelemental_water):
    obmol = cheminfo.qcelemental2OBMol(qcelemental_water)
    assert obmol.NumAtoms() == 3
