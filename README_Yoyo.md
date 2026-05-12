# CMB114 CW3 Molecule viewer and spectroscopy tool

## Installation
User is expected to install the [CIRpy library](https://cirpy.readthedocs.io/en/latest/) and [RDKit library](https://www.rdkit.org/docs/GettingStartedInPython.html) for drawing molecules. Internet is needed when accessing CIRpy.

User also expected to install [matplotlib](https://matplotlib.org/stable/index.html), but it usually installed when running Python.

## Mode - Molecule viewer
The viewer ables to generate the 2D and 3D molecule structure.
When entering a molecule that is from database, its info shows, such as molecular weight, spectra peaks and more. The strucutre generates from RDkit. For molecules outside of the database, no info is able to show. Their structures require to get SMILES string from CIRpy and therefore uses for generate the structures, which required installation and internet. If internet or CIRpy is unavailable, an error message will pop up.

View 2D structure:
1. Enter the chemical name.
2. Click "Draw" button.
3. The structure shows and export a PNG file of it.

View 3D structure:
1. Enter the chemical name.
2. Click "Draw" button.
3. Choose to export Mol or XYZ file.
4. Open the file in Avogadro to view the 3D sturcture.

## Mode - Molecule comparator
The comparator ables to compare the spectra of two of the molecules among 15 molecules from the database.

1. Enter the two molecules you want to compare and click "Compare".
2. IR, NMR, MS, and UV spectra of both molecules are shown. In summary box, a summary and comparison from the database of the molecules are also shown, such as their formula, molecular weight, and functional group.

## References for the Modes
Plotting with matplotlib: https://matplotlib.org/stable/plot_types/index.html
Grid layout of the subplots by matplotlib: https://matplotlib.org/stable/api/_as_gen/matplotlib.gridspec.GridSpec.html 
