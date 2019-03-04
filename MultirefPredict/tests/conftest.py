"""
Fixtures shared by tests in this directory
"""
import pytest
import qcelemental

@pytest.fixture(scope="module")
def qcelemental_water():
    mol = qcelemental.models.Molecule.from_data("""
    O                 0.000000000000     0.000000000000    -0.068516245955
    H                 0.000000000000    -0.790689888800     0.543701278274
    H                 0.000000000000     0.790689888800     0.543701278274
    """)
    return mol

@pytest.fixture(scope="module")
def xyz_water():
    xyz="""3

    O                 0.000000000000     0.000000000000    -0.068516245955
    H                 0.000000000000    -0.790689888800     0.543701278274
    H                 0.000000000000     0.790689888800     0.543701278274\n"""
    return xyz
