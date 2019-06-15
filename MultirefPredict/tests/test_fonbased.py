"""
Unit and regression test for the MultirefPredict package.
"""

# Import package, test suite, and other packages as needed
import pytest
from MultirefPredict.fonbased_diagnostic import FonBased
from qcengine.testing import using_terachem
from qcelemental.testing import compare_recursive
import pickle
from contextlib import contextmanager

@pytest.fixture(scope="class")
def fon_cu_complex(qcelemental_cu_complex):
    fonbased = FonBased(molecule=qcelemental_cu_complex, program="terachem")
    return fonbased

@pytest.fixture(scope="class")
def fon_trityl_radical(qcelemental_trityl_radical):
    fonbased = FonBased(molecule=qcelemental_trityl_radical, program="terachem")
    return fonbased

@pytest.fixture(scope="class")
def trityl_radical_task(datadir):
    task_file_path = datadir+"fon_trityl_task.obj"
    with open(task_file_path, 'rb') as task_file:
        task = pickle.load(task_file)
    return task

@pytest.fixture(scope="class")
def fon_trityl_result(datadir):
    file_path = datadir + "fon_trityl_result.obj"
    with open(file_path, 'rb') as result_file:
        result = pickle.load(result_file)
    return result

@pytest.fixture(scope="class")
def fon_cu_result(datadir):
    file_path = datadir + "fon_cu_result.obj"
    with open(file_path, 'rb') as result_file:
        result = pickle.load(result_file)
    return result

def test_fon_init(qcelemental_trityl_radical):
    fonbased = FonBased(molecule=qcelemental_trityl_radical, program="terachem")
    assert fonbased.fons is None
    assert fonbased.restricted == False
    assert fonbased.norb == 0
    assert fonbased.ncore == 0
    assert fonbased.nele == 0

@pytest.mark.parametrize("xc_list, expected_temperature", [
    ("pbe", 0.0158),
    ("blyp", 0.0158),
    ("svwn", 0.0158),
    ("revpbe", 0.0158),
    ("bop", 0.0158),
    ("b97", 0.0158),
    ("b3lyp", 0.0285),
    ("b3lyp1", 0.0285),
    ("b3lyp3", 0.0285),
    ("b3lyp5", 0.0285),
    ("bhandhlyp", 0.0475),
    ("pbe0", 0.0316),
    ("revpbe0", 0.0316),
    ])
def test_fon_setTemperature(xc_list, expected_temperature, fon_cu_complex):
    Thre = 1e-3
    assert compare_recursive(fon_cu_complex.setTemperature(xc_list), expected_temperature, atol = Thre)

@contextmanager
def does_not_raise():
    yield

@pytest.mark.parametrize("xc, expectation", [
    ("pbe", does_not_raise()),
    ("blyp", does_not_raise()),
    ("svwn", does_not_raise()),
    ("revpbe", does_not_raise()),
    ("bop", does_not_raise()),
    ("b97", does_not_raise()),
    ("b3lyp", does_not_raise()),
    ("b3lyp1", does_not_raise()),
    ("b3lyp3", does_not_raise()),
    ("b3lyp5", does_not_raise()),
    ("bhandhlyp", does_not_raise()),
    ("pbe0", does_not_raise()),
    ("revpbe0", does_not_raise()),
    ("wpbeh", pytest.raises(ValueError)),
    ("wpbe", pytest.raises(ValueError)),
    ("camb3lyp", pytest.raises(ValueError)),
    ])
def test_fon_setTemperature_exception(xc, expectation, fon_cu_complex):
    with expectation:
        fon_cu_complex.setTemperature(xc)

@using_terachem
def test_fon_computeFon(fon_cu_complex):
    result = fon_cu_complex.computeFon("PBE")
    assert fon_cu_complex.fons is not None

def test_fon_task(fon_trityl_radical, trityl_radical_task):
    with pytest.raises(ValueError):
        fon_task = fon_trityl_radical.FonTask(fon_trityl_radical.molecule, "psi4", "PBE")

    fon_task = fon_trityl_radical.FonTask(fon_trityl_radical.molecule,"terachem","PBE")
    Thre = 1e-6
    assert compare_recursive(fon_task.dict()['keywords'], trityl_radical_task.dict()['keywords'], atol=Thre)

def test_harvestFon(fon_trityl_radical, fon_trityl_result, fon_cu_complex, fon_cu_result):
    #first test trityl_radical
    fon_trityl_radical.FonTask(fon_trityl_radical.molecule, fon_trityl_radical.program, "PBE")
    fons=fon_trityl_radical.harvestFon(fon_trityl_result)
    expected_nele = 129
    assert fon_trityl_radical.nele == expected_nele
    assert fon_trityl_radical.restricted == False
    assert len(fons) > 0
    #Then test cu_complex
    fon_cu_complex.FonTask(fon_cu_complex.molecule, fon_cu_complex.program, "PBE")
    fons = fon_cu_complex.harvestFon(fon_cu_result)
    expected_nele = 120
    assert fon_cu_complex.nele == expected_nele
    assert fon_cu_complex.restricted == True
    assert len(fons) > 0

@using_terachem
def test_fon_computeDiagnostic(fon_cu_complex):
    Thre = 1e-2
    diag = fon_cu_complex.computeDiagnostic()
    expected = {'FOD': 0.980722,
                'Mattito':{ 'I_D': 0.473731315028, 'I_ND':  0.431365088635}, 
                'Entanglement': 0.732773215873}
    print(diag)
    assert compare_recursive(diag, expected, atol=Thre)

@using_terachem
def test_fon_computeFon_unrestricted(fon_trityl_radical):
    molecule_task = fon_trityl_radical.FonTask(fon_trityl_radical.molecule,\
                    fon_trityl_radical.program, "PBE")
    assert molecule_task is not None
    result = fon_trityl_radical.computeFon("PBE")
    assert fon_trityl_radical.fons is not None
    with pytest.raises(ValueError):
        failed_result = fon_trityl_radical.computeFon("TPSS")

@using_terachem
def test_fon_unrestricted(fon_trityl_radical):
    Thre = 1e-2
    diag = fon_trityl_radical.computeDiagnostic()
    expected = {'FOD': 0.34399500000000005, 
                'Mattito': {'I_D': 0.39268800205639215, 'I_ND': 0.16333423895850002},
                'Entanglement': -1}

    print(diag)
    assert compare_recursive(diag, expected, atol=Thre)

@using_terachem
def test_fon_b3lyp(fon_trityl_radical):
    Thre = 1e-2
    diag2 = fon_trityl_radical.computeDiagnostic("b3lyp")
    expected2 = {'FOD': 0.7340689999999996, 
                 'Mattito': {'I_D': 0.7479020997602953, 'I_ND': 0.33741997473349994},
                 'Entanglement': -1}

    print(diag2)
    assert compare_recursive(diag2, expected2, atol=Thre)
