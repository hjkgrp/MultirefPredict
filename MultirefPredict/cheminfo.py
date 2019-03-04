import qcelemental
import openbabel

def qcelemental2xyz(molecule):
    coords = molecule.pretty_print().split('\n')[4:-2]
    xyz = str(len(coords)) + "\n\n"
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
    mol_str = str(charge)+" "+str(spinmult)
    for i in range(2, len(xyz)):
        mol_str += xyz[i]
    molecule = qcelemental.models.Molecule.from_data(mol_str)
    return molecule

