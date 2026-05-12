# Code by Mo
import os

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

#  basically looks for peak clusters within the range of the groups of the molecule you typed in
def infer_functional_groups(ir_peaks):
    #checking for clusters of peaks
    def _has(lo, hi, min_i=0):
        for wn, inten, label in ir_peaks:
            if lo <= wn <= hi and inten >= min_i:
                return True
        return False

    inferences = []
    reasoning  = []

# alcohol
    if _has(3200, 3550) and _has(900, 1300):
        if   _has(1020, 1075): inferences.append("primary alcohol");   reasoning.append("O-H + C-O ~1050 -> primary alcohol")
        elif _has(1100, 1150): inferences.append("secondary alcohol"); reasoning.append("O-H + C-O ~1125 -> secondary alcohol")
        elif _has(1150, 1210): inferences.append("tertiary alcohol");  reasoning.append("O-H + C-O ~1175 -> tertiary alcohol")
        else:                  inferences.append("alcohol");           reasoning.append("Broad O-H -> alcohol or phenol")

# carboxylic acid
    if _has(2500, 3300) and _has(1700, 1725):
        inferences.append("carboxylic acid"); reasoning.append("Broad O-H (2500-3300) + C=O ~1710 -> carboxylic acid")

# primary  and secondary amine
    if _has(3350, 3500) and _has(3280, 3380):
        inferences.append("primary amine"); reasoning.append("Two N-H bands -> primary amine")
    elif _has(3300, 3450) and not _has(3200, 3550):
        inferences.append("secondary amine"); reasoning.append("Single N-H -> secondary amine")

# aldehyde
    if _has(2700, 2750) and _has(2800, 2870):
        inferences.append("aldehyde"); reasoning.append("Aldehyde C-H doublet (~2720, ~2820) -> aldehyde")

#alkyne
    if _has(2100, 2260): inferences.append("alkyne"); reasoning.append("C≡C (2100-2260) -> alkyne")
    if _has(2100, 2200) and "alkyne" not in inferences:
        inferences.append("nitrile"); reasoning.append("C≡N -> nitrile")
# different carbonyl compounds
    if   _has(1800, 1850): inferences.append("acyl halide"); reasoning.append("C=O (1800-1850) -> acyl halide")
    elif _has(1730, 1800): inferences.append("ester");       reasoning.append("C=O (1730-1800) -> ester")
    elif _has(1700, 1730) and "carboxylic acid" not in inferences and "aldehyde" not in inferences:
        inferences.append("ketone"); reasoning.append("C=O (1700-1730) -> ketone")
    elif _has(1630, 1700):
        if "amine" in " ".join(inferences): inferences.append("amide");  reasoning.append("C=O + N-H -> amide")
        elif _has(1600, 1660):              inferences.append("alkene"); reasoning.append("C=C (1600-1680) -> alkene")

#aromatic compounds
    if _has(1580, 1620) and _has(1470, 1510):
        inferences.append("aromatic"); reasoning.append("C=C aromatic (~1600, ~1500) -> aromatic ring")
        if _has(690, 750): inferences.append("monosubstituted benzene")

    if _has(1100, 1300, 70) and not any(k in inferences for k in ["alcohol","ester","ketone","aldehyde"]):
        inferences.append("ether (possible)")

    if _has(600, 800, 60): inferences.append("alkyl chloride (possible)")
    if _has(500, 600, 50): inferences.append("alkyl bromide (possible)")

# to avoid any duplicates or repetitions
    seen = set()
    unique = []
    for x in inferences:
        if x not in seen:
            seen.add(x)
            unique.append(x)
    return unique, reasoning

# code written by Yoyo
def draw_molecule_2d(molecule, save_path=None, size=(400, 300)):
    # generates a 2D structure PNG using RDKit
    rdkit_mol = molecule.get_rdkit_mol()
    if rdkit_mol is None:
        return None

    AllChem.Compute2DCoords(rdkit_mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    drawer.drawOptions().addAtomIndices = False
    drawer.drawOptions().addStereoAnnotation = True
    drawer.DrawMolecule(rdkit_mol)
    drawer.FinishDrawing()

    if save_path is None:
        os.makedirs("output", exist_ok=True)
        save_path = f"output/{molecule.name.replace(' ', '_')}_2D.png"

    with open(save_path, "wb") as f:
        f.write(drawer.GetDrawingText())

    return save_path
# end of code written Yoyo