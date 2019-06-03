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

@pytest.fixture(scope="class")
def fon_trityl_radical(qcelemental_trityl_radical):
    fonbased = FonBased(molecule=qcelemental_trityl_radical, program="terachem")
    return fonbased
@pytest.fixture(scope="class")
def trityl_radical_task():
    task={
     'molecule': {'symbols': ['C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'C',    'H',    'H',    'H',    'H',    'H',    'H',    'H',    'H',    'H',    'H',    'H',    'H',    'H',    'H',    'H'],
      'geometry': [-2.96111532,    0.28836487,    -0.76910104,    -1.94737287,    1.54462598,    1.33822583,    0.62061102,    1.33578573,    1.91269075,    2.27278214,    -0.15571369,    0.40314734,    1.20125467,    -1.41463107,    -1.71789112,    -1.3655295,    -1.1885923,    -2.29190052,    4.96173057,    -0.38595766,    1.01423875,    6.30474663,    1.74610788,    2.15876124,    5.66330434,    4.2883185,    1.56051828,    6.94786901,    6.31857978,    2.65492965,    8.91069273,    5.89210824,    4.39004512,    9.57424837,    3.39852053,    5.01308084,    8.30293209,    1.36144981,    3.91599397,    6.30409904,    -2.74828654,    0.49320078,    8.90714342,    -2.76239181,    -0.18491704,    10.17991229,    -5.01956511,    -0.68771999,    8.91237483,    -7.34595324,    -0.52311453,    6.34397486,    -7.37998961,    0.14779566,    5.05923769,    -5.12783951,    0.64237347,    -4.9733875,    0.45707704,    -1.21876948,    -3.17489669,    2.68321734,    2.55528504,    1.38329522,    2.29258396,    3.57944432,    2.43245727,    -2.54023328,    -2.9397217,    -2.1282345,    -2.16077959,    -3.95249926,    4.16393795,    4.64300787,    0.18176998,    6.42559692,    8.25332346,    2.136773,    9.91174512,    7.48595401,    5.2489155,    11.08263202,    3.0392172,    6.38405442,    8.80687823,    -0.57191479,    4.44857172,    9.90370775,    -0.9597762,    -0.36508343,    12.17596071,    -4.96444584,    -1.23307949,    9.91450146,    -9.113073,    -0.91418413,    5.3391323,    -9.18272773,    0.30504866,    3.07140552,    -5.18063471,    1.20859218],
      'masses': [12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    12.0,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223,    1.00782503223],
      'fragment_charges': [0.0],
      'fragment_multiplicities': [2],
      'schema_name': 'qcschema_molecule',
      'schema_version': 2,
      'name': 'C19H15',
      'molecular_charge': 0.0,
      'molecular_multiplicity': 2,
      'atomic_numbers': [6,    6,    6,    6,    6,    6,    6,    6,    6,    6,    6,    6,    6,    6,    6,    6,    6,    6,    6,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
      'mass_numbers': [12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    12,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1],
      'fix_com': False,
      'fix_orientation': False,
      'provenance': {'creator': 'QCElemental',
       'version': 'v0.4.0+6.ga9026b6',
       'routine': 'qcelemental.molparse.from_string'}},
     'model': {'method': 'uPBE', 'basis': 'lacvps_ecp'},
     'id': None,
     'schema_name': 'qcschema_input',
     'schema_version': 1,
     'keywords': {'gpus': '1',
      'maxit': '1500',
      'scf': 'diis+a',
      'convthre': '1e-6',
      'precision': 'double',
      'units': 'bohr',
      'fon': 'yes',
      'fon_method': 'fermi',
      'fon_temperature': 0.0158,
      'fon_print': '1',
      'method': 'uPBE',
      'closed': 19,
      'active': 296},
     'extras': {},
     'provenance': {'creator': 'QCElemental',
      'version': 'v0.4.0+6.ga9026b6',
      'routine': 'qcelemental.models.results'}}
    return task

def test_fon_init(qcelemental_trityl_radical):
    fonbased = FonBased(molecule=qcelemental_trityl_radical, program="terachem")
    assert fonbased.fons is None
    assert fonbased.restricted == False
    assert fonbased.norb == 0
    assert fonbased.ncore == 0
    assert fonbased.nele == 0

@using_terachem
def test_fon_computeFon(fon_cu_complex):
    result = fon_cu_complex.computeFon("PBE")
    assert fon_cu_complex.fons is not None

def test_fon_task(fon_trityl_radical, qcelemental_trityl_radical, trityl_radical_task):
    with pytest.raises(ValueError):
        fon_task = fon_trityl_radical.FonTask(qcelemental_trityl_radical, "psi4", "PBE")

    fon_task = fon_trityl_radical.FonTask(qcelemental_trityl_radical,"terachem","PBE")
    Thre = 1e-6
    assert compare_recursive(fon_task.dict()['keywords'], trityl_radical_task['keywords'], atol=Thre)


@using_terachem
def test_fon_computeDiagnostic(fon_cu_complex):
    Thre = 1e-3
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

@using_terachem
def test_fon_unrestricted(fon_trityl_radical):
    Thre = 1e-3
    diag = fon_trityl_radical.computeDiagnostic()
    expected = {'FOD': 0.34399500000000005, 
                'Mattito': {'I_D': 0.39268800205639215, 'I_ND': 0.16333423895850002},
                'Entanglement': -1}

    print(diag)
    assert compare_recursive(diag, expected, atol=Thre)
