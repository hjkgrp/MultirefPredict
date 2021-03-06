"""
Fixtures shared by tests in this directory
"""
import pytest
import qcelemental
from pkg_resources import resource_filename, Requirement

@pytest.fixture(scope="module")
def qcelemental_water():
    mol = qcelemental.models.Molecule.from_data("\n"\
       + "O                 0.000000000000     0.000000000000    -0.068516245955\n"\
       + "H                 0.000000000000    -0.790689888800     0.543701278274\n"\
       + "H                 0.000000000000     0.790689888800     0.543701278274\n"\
       +"\n"
    )
    return mol

@pytest.fixture(scope="module")
def xyz_water():
    xyz="3\n\n"\
       + "    O                 0.000000000000     0.000000000000    -0.068516245955\n"\
       + "    H                 0.000000000000    -0.790689888800     0.543701278274\n"\
       + "    H                 0.000000000000     0.790689888800     0.543701278274\n"
    return xyz

@pytest.fixture(scope="module")
def qcelemental_cu_complex():
    mol = qcelemental.models.Molecule.from_data("\n"\
        +"2 1\n"\
        +"O -0.420275 1.402468 0.101942\n"\
        +"O -0.440639 -0.486438 -0.102059\n"\
        +"Cu 1.098238 0.442692 -0.009679\n"\
        +"Cu -1.959153 0.473350 0.009562\n"\
        +"N -3.386852 1.843664 0.160103\n"\
        +"N -3.413724 -0.868627 -0.145405\n"\
        +"N 2.552792 1.784687 0.145285\n"\
        +"N 2.525952 -0.927604 -0.160217\n"\
        +"C 3.838802 1.128848 -0.261856\n"\
        +"C 3.818294 -0.300514 0.271356\n"\
        +"C -4.679199 1.216590 -0.271474\n"\
        +"C -4.699725 -0.212774 0.261735\n"\
        +"H 2.605864 2.113287 1.123426\n"\
        +"H 2.363065 2.624052 -0.423452\n"\
        +"H 2.587660 -1.250217 -1.139828\n"\
        +"H 2.310711 -1.766880 0.399591\n"\
        +"H 4.714722 1.690265 0.115054\n"\
        +"H 3.882273 1.135693 -1.367752\n"\
        +"H 3.840670 -0.308405 1.377933\n"\
        +"H 4.688635 -0.881164 -0.089123\n"\
        +"H -3.171601 2.682938 -0.399701\n"\
        +"H -3.448557 2.166273 1.139714\n"\
        +"H -3.224006 -1.707993 0.423331\n"\
        +"H -3.466797 -1.197225 -1.123545\n"\
        +"H -4.701574 1.224483 -1.378049\n"\
        +"H -5.549535 1.797248 0.089003\n"\
        +"H -5.575651 -0.774178 -0.115176\n"\
        +"H -4.743200 -0.219619 1.367630\n"\
        +"\n"
    )
    return mol

@pytest.fixture(scope="module")
def qcelemental_trityl_radical():
    mol = qcelemental.models.Molecule.from_data("\n"\
        + "0 2\n"\
        + "C    -1.5669547451    0.1525961178   -0.4069907437\n"\
        + "C    -1.0305053421    0.8173808668    0.7081586097\n"\
        + "C     0.3284132085    0.7068673641    1.0121523585\n"\
        + "C     1.2027045139   -0.0824001370    0.2133363853\n"\
        + "C     0.6356765971   -0.7485905261   -0.9090688316\n"\
        + "C    -0.7226070902   -0.6289759573   -1.2128215257\n"\
        + "C     2.6256347453   -0.2042399959    0.5367120349\n"\
        + "C     3.3363282370    0.9240004987    1.1423672507\n"\
        + "C     2.9968915962    2.2692804227    0.8257907091\n"\
        + "C     3.6766539433    3.3436484234    1.4049282644\n"\
        + "C     4.7153355259    3.1179694060    2.3231118289\n"\
        + "C     5.0664740448    1.7984196135    2.6528081376\n"\
        + "C     4.3937224416    0.7204482120    2.0722547643\n"\
        + "C     3.3359855438   -1.4543306062    0.2609906110\n"\
        + "C     4.7134573125   -1.4617947915   -0.0978538809\n"\
        + "C     5.3869775886   -2.6562394624   -0.3639257485\n"\
        + "C     4.7162256544   -3.8873110456   -0.2768202853\n"\
        + "C     3.3570869204   -3.9053223190    0.0782100962\n"\
        + "C     2.6772332870   -2.7135358088    0.3399294001\n"\
        + "H    -2.6318033263    0.2418747524   -0.6449450360\n"\
        + "H    -1.6800829740    1.4198974700    1.3521986109\n"\
        + "H     0.7320083085    1.2131831851    1.8941603601\n"\
        + "H     1.2872009507   -1.3442335612   -1.5556337305\n"\
        + "H    -1.1262131968   -1.1434353152   -2.0915725316\n"\
        + "H     2.2034610687    2.4569739529    0.0961885326\n"\
        + "H     3.4002794568    4.3674706849    1.1307315735\n"\
        + "H     5.2450696331    3.9613962639    2.7776064628\n"\
        + "H     5.8646763007    1.6082844786    3.3782961091\n"\
        + "H     4.6603992579   -0.3026442721    2.3540827763\n"\
        + "H     5.2408164435   -0.5078916932   -0.1931938322\n"\
        + "H     6.4432409235   -2.6270716002   -0.6525175643\n"\
        + "H     5.2465282301   -4.8224305525   -0.4837654065\n"\
        + "H     2.8253471361   -4.8592902457    0.1614248002\n"\
        + "H     1.6253178047   -2.7414738227    0.6395594407\n"\
        + "\n"
    )
    return mol
@pytest.fixture(scope="module")
def datadir():
    mydir = resource_filename(Requirement.parse("MultirefPredict"),"MultirefPredict/tests/data/")
    return mydir
