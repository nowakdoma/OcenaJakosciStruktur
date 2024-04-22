# Structure Quality Assessment

This repository contains two Python scripts for assessing the quality of protein structures: `superpozycja.py` and `clashscore.py`. These scripts utilize the Biopython package for parsing PDB files and performing structural analyses.

## RMSD

### Description:
This script calculates the Root Mean Square Deviation (RMSD) between a reference structure and a model structure. It superimposes the model structure onto the reference structure to minimize the RMSD.

### Features:
- **RMSD Calculation:** Calculates the RMSD between the reference and model structures.
- **Superimposition:** Superimposes the model structure onto the reference structure for RMSD calculation.
- **Identification of Best and Worst Models:** Identifies the model with the lowest and highest RMSD values.

### Usage:
1. Ensure Python and the required libraries are installed.
2. Run the script from the command line.
3. Provide the path to the reference structure using the `--r` argument.
4. Provide the path to the model structure using the `--m` argument.

#### Example:
```
python superpozycja.py --r reference.pdb --m model.pdb
```
### Output:
- The script prints the RMSD values for each model.
- It saves the model with the lowest RMSD as `best.pdb` and the model with the highest RMSD as `worst.pdb`.

## Clashscore

### Description:
This script calculates the Clash Score for a protein structure. The Clash Score indicates the number of atomic clashes or overlaps within the structure.

### Features:
- **Clash Score Calculation:** Calculates the Clash Score based on atomic clashes within the structure.

### Usage:
1. Ensure Python and the required libraries are installed.
2. Run the script from the command line.
3. Provide the path to the PDB file containing the protein structure using the `--f` argument.

#### Example:
```
python clashscore.py --f protein_structure.pdb
```
### Output:
- The script prints the calculated Clash Score for the provided protein structure.

## Dependencies:
- Python 3.x
- Biopython

## Note:
- Ensure proper installation of Biopython before running the scripts.
- These scripts assume that the input PDB files follow the standard PDB format.
