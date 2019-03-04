import qcelemental
import openbabel

def qcelemental2Xyz(molecule):
    coords = molecule.pretty_print().split('\n')[4:-2]
    xyz = str(len(coords)) + "\n\n"
    for line in coords:
        xyz += line + "\n"
    return xyz

def qcelemental2OBMol(molecule):
    xyz = qcelemental2Xyz(molecule)
    obmol = openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz","smi")
    obConversion.ReadString(obmol,xyz)
    return obmol
