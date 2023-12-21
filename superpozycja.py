import argparse
from Bio import PDB
from Bio.PDB import Superimposer
import numpy as np

def parseArguments():
    parser = argparse.ArgumentParser(description='Calculate RMSD between reference and model structures.')
    parser.add_argument('--r', required=True, help='Path to the reference file')
    parser.add_argument('--m', required=True, help='Path to the model file')
    return parser.parse_args()

def getStructure(filePath, structureId):
    parser = PDB.PDBParser(QUIET=True)
    return parser.get_structure(structureId, filePath)

def getMatchingAtoms(referenceAtoms, modelAtoms):
    refAtoms = []
    subAtoms = []

    for refAtom in referenceAtoms:
        for subAtom in modelAtoms:
            if refAtom.get_parent().get_id()[1] == subAtom.get_parent().get_id()[1] and refAtom.get_id() == subAtom.get_id():
                refAtoms.append(refAtom)
                subAtoms.append(subAtom)

    return refAtoms, subAtoms

def superimposeAtoms(referenceAtoms, modelAtoms):
    superimposer = Superimposer()
    superimposer.set_atoms(referenceAtoms, modelAtoms)
    superimposer.apply(modelAtoms)

def calculateRmsd(referenceAtoms, modelAtoms):
    modelCoord = np.array([atom.coord for atom in modelAtoms])
    referenceCoord = np.array([atom.coord for atom in referenceAtoms])
    diff = modelCoord - referenceCoord
    sumOfSquares = np.sum(diff ** 2)
    return np.sqrt(sumOfSquares / len(referenceCoord))

def writePDB(filename, structure):
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(filename)

def main():
    args = parseArguments()

    referenceStructure = getStructure(args.r, 'reference')
    modelStructure = getStructure(args.m, 'structures')

    models = [model for model in modelStructure]

    rmsdValues = []

    for model in models:
        referenceAtoms = [atom for atom in referenceStructure.get_atoms()]
        modelAtoms = [atom for atom in model.get_atoms()]

        refAtoms, subAtoms = getMatchingAtoms(referenceAtoms, modelAtoms)
        superimposeAtoms(refAtoms, subAtoms)

        rmsd = calculateRmsd(refAtoms, subAtoms)
        rmsdValues.append(rmsd)
        print(rmsd)


    bestModel = models[np.argmin(rmsdValues)]
    worstModel = models[np.argmax(rmsdValues)]

    writePDB('best.pdb', bestModel)
    writePDB('worst.pdb', worstModel)

if __name__ == "__main__":
    main()
