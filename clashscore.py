import argparse
from Bio.PDB import PDBParser, NeighborSearch

def parseArguments():
    parser = argparse.ArgumentParser(description='Calculate Clash Score for a PDB structure.')
    parser.add_argument('--f', required=True, help='Path to the PDB file')
    return parser.parse_args()

def calculateClashScore(structure):
    vanDerWaalsRadii = {"H": 1.2, "O": 1.52, "N": 1.55, "C": 1.7}
    atoms = list(structure.get_atoms())
    numAtoms = len(atoms)
    visitedPairs = set()
    neighborSearch = NeighborSearch(atoms)
    badOverlaps = 0

    for atom1 in atoms:
        res1 = atom1.get_parent()
        for atom2 in neighborSearch.search(atom1.coord, 4, level='A'):
            res2 = atom2.get_parent()

            atomPair = frozenset([atom1, atom2])

            if atomPair in visitedPairs:
                continue
            visitedPairs.add(atomPair)

            if abs(res1.id[1] - res2.id[1]) < 2 or res1 == res2:
                continue

            distance = atom1 - atom2

            radius1 = vanDerWaalsRadii.get(atom1.element, 1.5)
            radius2 = vanDerWaalsRadii.get(atom2.element, 1.5)

            if (radius1 + radius2 - distance) >= 0.4:
                badOverlaps += 1

    clashScore = (1000 * badOverlaps) / numAtoms
    return clashScore

def main():
    args = parseArguments()
    pdbParser = PDBParser(QUIET=True)
    structure = pdbParser.get_structure("model", args.f)
    clashScore = calculateClashScore(structure)
    print("Clash Score:", clashScore)

if __name__ == "__main__":
    main()
