# By Yoyo
import json
import os
from dataclasses import dataclass, field

# try to import CIRpy library for drawing the molecule
try:
    import cirpy
    CIRPY_AVAILABLE = True
except ImportError:
    CIRPY_AVAILABLE = False


@dataclass
class Molecule:
    """
    define data type
    """
    name: str
    iupac_name: str
    formula: str
    molecular_weight: float
    smiles: str
    ir_peaks: list = field(default_factory=list)
    ms_peaks: list = field(default_factory=list)
    uv_peaks:  list = field(default_factory=list)
    nmr_peaks: list = field(default_factory=list)
    functional_groups: list = field(default_factory=list)
    compound_class: str = ""
    fun_fact: str = ""

    def get_rdkit_mol(self):
        """
        construct molecule from a SMILES string
        """
        from rdkit import Chem
        mol = Chem.MolFromSmiles(self.smiles) 
        if mol is None:
            print(f"[WARNING] Bad SMILES for {self.name}")
        return mol

    def get_base_peak(self):
        return max(self.ms_peaks, key=lambda p: p[1]) if self.ms_peaks else None

    def get_molecular_ion(self):
        return next((mz for mz, _, __ in self.ms_peaks if abs(mz - self.molecular_weight) <= 1.0), None)

    def summary(self):
        return (f"{self.name.title()} ({self.formula}, MW={self.molecular_weight} g/mol)\n"
                f"  IUPAC: {self.iupac_name}\n  Class: {self.compound_class}\n"
                f"  Groups: {', '.join(self.functional_groups) or 'none'}")


class MoleculeDatabase:
    def __init__(self, json_path=None):
        self._compounds = {}
        if json_path is None:
            json_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "molecules.json")
        self._load_from_json(json_path)
        print(f"[Database] {len(self._compounds)} compounds loaded.")
    
    def _load_from_json(self, path):
        """
        import data from database in JSON file
        """
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        for entry in data["molecules"]:

            # showing the info of the molecule
            mol = Molecule(
                name = entry["name"],
                iupac_name = entry["iupac_name"],
                formula = entry["formula"],
                molecular_weight  = entry["molecular_weight"],
                smiles = entry["smiles"],
                compound_class = entry.get("compound_class", ""),
                functional_groups = entry.get("functional_groups", []),
                fun_fact  = entry.get("fun_fact", ""),
                ir_peaks = [tuple(p) for p in entry.get("ir_peaks",  [])],
                ms_peaks = [tuple(p) for p in entry.get("ms_peaks",  [])],
                uv_peaks = [tuple(p) for p in entry.get("uv_peaks",  [])],
                nmr_peaks = [tuple(p) for p in entry.get("nmr_peaks", [])],
            )
            self._compounds[mol.name.lower()] = mol

    def get(self, name):
        return self._compounds.get(name.lower())

    def get_or_fetch(self, name):
        mol = self.get(name)
        if mol:
            return mol, True
        if CIRPY_AVAILABLE:
            try:
                # get SMILES string of the molecule from CIRpy
                smiles = cirpy.resolve(name, "smiles") 
                
                # no data available for the molecule generated with CIRpy
                if smiles:
                    mol = Molecule(
                        name=name.lower(), iupac_name=name, formula="Unknown", molecular_weight=0.0, smiles=smiles,                       
                        ir_peaks=[], ms_peaks=[], uv_peaks=[], nmr_peaks=[], functional_groups=[],
                        compound_class="Unknown", fun_fact=""
                    )
                    return mol, False
            except Exception:
                pass
        return None, False

    def all_molecules(self):
        return list(self._compounds.values())

    def list_names(self):
        return sorted(self._compounds.keys())
    
