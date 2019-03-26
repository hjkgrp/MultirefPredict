import qcelemental
import openbabel

def qcelemental2xyz(molecule):
    mol_str = molecule.pretty_print().strip('\n').split('\n')
    natom = len(molecule.symbols)
    startline = -natom
    coords = mol_str[startline:]
    xyz = str(natom) + "\n\n"
    for line in coords:
        xyz += line + "\n"
    return xyz

def qcelemental2OBMol(molecule):
    xyz = qcelemental2xyz(molecule)
    obmol = openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz","smi")
    obConversion.ReadString(obmol,xyz)
    return obmol

def xyzfile2qcelemental(xyzfile, charge=0, spinmult=1):
    xyz = open(xyzfile).readlines()
    mol_str = str(charge)+" "+str(spinmult) + "\n"
    for i in range(2, len(xyz)):
        mol_str += xyz[i]
    molecule = qcelemental.models.Molecule.from_data(mol_str)
    return molecule
