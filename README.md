# CMB114 CW3 Molecule viewer and spectroscopy tool by Yoyo and Muhammed

Contributions of code are defined in each file.

You can run the code by calling:

~~~~
python main_gui.py
~~~~

## General info
This tool aims to view the molecule structure and its spectra data accessed in the database. User can view any molecule structure when CIRpy is available, but without its spectra data.

## Installation
User is expected to install the [CIRpy library](https://cirpy.readthedocs.io/en/latest/) and [RDKit library](https://www.rdkit.org/docs/GettingStartedInPython.html) for drawing molecules. Internet is needed when accessing CIRpy.

User also expected to install [matplotlib](https://matplotlib.org/stable/index.html), but it usually installed when running Python.

## Mode - Predictor
Type a molecule name from the database (ethanol, acetone, etc) and select its respective spectra to generate, either IR, NMR, Mass-spec, UV-VIS or all four. Data was sourced using the sdbs website.
Using figurecanvas from tkinter made embedding the graphs onto the UI more straightforward.
It is worth noting that the user's scroll wheel may not register as this was an issue that couldnt be addressed in time, a scrollable canvas and scroll bar was used in the mean time.

## Sub-Mode - Database
The viewer gets to see the compounds already built in and a fun fact about them. The option to click each one and redirect you to the predictor to generate the plots was suggested but couldn't be implemented in time. Worth looking into in the future.
The database originally started off as python code but after meetings with our project supervisor it was suggested we use JSON, so the work around was to use pythons json module to export to JSON file and use an online JSON formatter to get it in the right format and ensure it was less messy.

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
1. Plotting with matplotlib: https://matplotlib.org/stable/plot_types/index.html
2. Grid layout of the subplots by matplotlib: https://matplotlib.org/stable/api/_as_gen/matplotlib.gridspec.GridSpec.html
3. SDBS Database: https://sdbs.db.aist.go.jp
4. Scrollbar: https://tkdocs.com/tutorial/canvas.html
5. Embedding matplotlib into tkinter: https://matplotlib.org/3.1.0/gallery/user_interfaces/embedding_in_tk_sgskip.html
