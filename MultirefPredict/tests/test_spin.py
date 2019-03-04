"""
Unit and regression test for the MultirefPredict package.
"""

# Import package, test suite, and other packages as needed
import pytest
import sys
from MultirefPredict.spin import atomSpinMultDict

@pytest.mark.parametrize("symbol_list, expected_spinmult" , [
        ("H", 2),
        ("He", 1),
        ("Li", 2),
        ("Be", 1),
        ("B", 2),
        ("C", 3),
        ("N", 4),
        ("O", 3),
        ("F", 2),
        ("Ne", 1),
        ("Na", 2),
        ("Mg", 1),
        ("Al", 2),
        ("Si", 3),
        ("P", 4),
        ("S", 3),
        ("Cl", 2),
        ("Ar", 1),
        ("K", 2),
        ("Ca", 1),
        ("Sc", 2),
        ("Ti", 3),
        ("V", 4),
        ("Cr", 7),
        ("Mn", 6),
        ("Fe", 5),
        ("Co", 4),
        ("Ni", 3),
        ("Cu", 2),
        ("Zn", 1),
        ("Ga", 2),
        ("Ge", 3),
        ("As", 4),
        ("Se", 3),
        ("Br", 2),
        ("Kr", 1),
])
def test_atomSpinMultDict(symbol_list, expected_spinmult):
    spindict = atomSpinMultDict()
    assert "MultirefPredict" in sys.modules
    assert spindict.get_spinmult(symbol_list) == expected_spinmult

@pytest.mark.parametrize("symbol_charge, spinmult" , [
    (["Sc", 2], 2),
    (["Ti", 2], 3),
    (["V", 2], 4),
    (["Cr", 2], 5),
    (["Mn", 2], 6),
    (["Fe", 2], 5),
    (["Co", 2], 4),
    (["Ni", 2], 3),
    (["Cu", 2], 2),
    (["Zn", 2], 1),
])
def test_atomSpinMultCation(symbol_charge,spinmult):
    spindict = atomSpinMultDict()
    assert spindict.get_spinmult(symbol_charge[0],symbol_charge[1]) == spinmult

