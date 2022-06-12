# RDKit tools
This repository containins tools built using RDKit to make life easier for cheminformatics.

## Scripts and their use:
1. **add-column-to-sdf-files-in-a-directory.py**: Can add a property to all SDF files in a folder with it's respective file name for easy backtracing.   Requires RDKit.
2. **excel-to-sdf-conversion.py**: Can convert an excel (.xlsx) file into an SDF file. Requires RDKit.
3. **random-sdf-selector.py**: Can select a given number of molecules from an SDF file randomly for unbiased sampling. Requires RDKit.
4. **synthetic_accessibility_scorer_SYBA.py**: Can calculate the synthetic accessibility scores for an SDF file. Requires RDKit and SYBA.
5. **E3FP-calculator.py**: Can calculate the E3FP fingerprints for molecules in an SDF file. Requires RDKit and E3FP.
6. **plip-tools**: Can generate pymol sessions for protein-ligand complexes in a folder (.pdb and .mae) formats.
