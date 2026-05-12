"""
main_gui.py — run with: python main_gui.py
"""

import tkinter as tk
from tkinter import ttk, messagebox
import os
from PIL import Image, ImageTk
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from molecules import MoleculeDatabase
from spectra import IRSpectrum, MassSpectrum, UVVisSpectrum, NMRSpectrum
from spectro_tools import draw_molecule_2d
from comparator import SpectrumComparator

# colours
BG = "#0F0F1A"; PANEL = "#1A1A2E"; CARD = "#16213E"
ACCENT = "#00BFFF"
TEXT = "#CCCCCC"; WHITE = "#FFFFFF"; GREEN = "#2ECC71"; GOLD = "#FFD700"

os.makedirs("output/spectra",    exist_ok=True)
os.makedirs("output/structures", exist_ok=True)

HELP_TEXT = """HOW TO USE THIS TOOL

PREDICTOR
  Enter a molecule name from the database.
  Select which spectra to generate (IR, MS, UV, NMR or all).
  The 2D structure appears, plots appear below it.
  Only works for the 15 built-in database molecules.

VIEWER
  Draw any molecule. Database molecules show full info.
  Unknown molecules are looked up online (needs internet).

COMPARATOR
  Enter two molecule names to compare side by side.
  Spectra are only available for database molecules.

DATABASE
  Browse all 15 built-in compounds.
"""

# code written by Yoyo and Mo
class SpectroscopyApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Spectroscopy Analysis Tool  CMB114 CW3")
        self.root.configure(bg=BG)
        self.root.geometry("1050x750")
        self.root.resizable(True, True)

        self.db = MoleculeDatabase()
        self.comparator = SpectrumComparator(self.db)
        self.viewer_mol  = None
        self.viewer_name = ""

        self._build_ui()
        self.show_predictor()

    def _build_ui(self):
        # header
        hdr = tk.Frame(self.root, bg=PANEL, pady=10)
        hdr.pack(fill="x")
        tk.Frame(hdr, bg=ACCENT, width=5).pack(side="left", fill="y")
        tk.Label(hdr, text="Spectroscopy Analysis Tool", bg=PANEL, fg=WHITE,
                 font=("Arial", 15, "bold")).pack(side="left", padx=14)
        tk.Label(hdr, text="CMB114 CW3  By Yoyo and Muhammed", bg=PANEL, fg=TEXT,
                 font=("Arial", 9)).pack(side="left", padx=4)
        self.make_btn(hdr, "?", self.show_help, side="right")

        # nav
        nav = tk.Frame(self.root, bg=BG, pady=8)
        nav.pack(fill="x")
        self.nav_btns = {}
        for label, cmd in [("Predictor", self.show_predictor),
                           ("Comparator", self.show_comparator),
                           ("Viewer", self.show_viewer),
                           ("Database", self.show_database)]:
            b = tk.Label(nav, text=label, bg=CARD, fg=WHITE,
                         font=("Arial", 10, "bold"), pady=8, cursor="hand2")
            b.pack(side="left", fill="x", expand=True, padx=2)
            b.bind("<Button-1>", lambda e, c=cmd: c())
            self.nav_btns[label] = b
        tk.Frame(self.root, bg=ACCENT, height=2).pack(fill="x")

        self.content = tk.Frame(self.root, bg=BG)
        self.content.pack(fill="both", expand=True)

    def make_btn(self, parent, text, cmd, side=None):
        # uses tk.Label with click binding so buttons never flash white
        b = tk.Label(parent, text=text, bg=ACCENT, fg=WHITE,
                     font=("Arial", 10, "bold"), pady=6, padx=12, cursor="hand2")
        b.bind("<Button-1>", lambda e: cmd())
        b.bind("<Enter>", lambda e: b.config(fg="#AAAAAA"))
        b.bind("<Leave>", lambda e: b.config(fg=WHITE))
        if side:
            b.pack(side=side, padx=(0, 6))
        else:
            b.pack(anchor="w", pady=6)
        return b

    def switch_to(self, label):
        for k, b in self.nav_btns.items():
            b.config(bg=ACCENT if k == label else CARD, fg=WHITE)
        self.root.update_idletasks()
        for widget in self.content.winfo_children():
            widget.destroy()
        plt.close("all")

    def scrollable(self, parent):
        # ref: tkdocs.com/tutorial/canvas.html
        outer = tk.Frame(parent, bg=BG)
        outer.pack(fill="both", expand=True)
        c = tk.Canvas(outer, bg=BG, highlightthickness=0)
        sb = ttk.Scrollbar(outer, orient="vertical", command=c.yview)
        inner = tk.Frame(c, bg=BG, padx=18, pady=10)
        inner.bind("<Configure>", lambda e: c.configure(scrollregion=c.bbox("all")))
        c.create_window((0, 0), window=inner, anchor="nw")
        c.configure(yscrollcommand=sb.set)
        c.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        def scroll(event):
            c.yview_scroll(int(-1 * (event.delta / 60)), "units")
        c.bind("<MouseWheel>", scroll)
        inner.bind("<MouseWheel>", scroll)
        return inner

    def show_help(self):
        pop = tk.Toplevel(self.root)
        pop.title("How to use")
        pop.configure(bg=PANEL)
        pop.geometry("480x320")
        pop.grab_set()
        tk.Label(pop, text=HELP_TEXT, bg=PANEL, fg=TEXT, font=("Courier", 9),
                 justify="left", wraplength=440).pack(padx=20, pady=14)
        self.make_btn(pop, "Close", pop.destroy)

    def embed_fig(self, p, fig):
        # embeds matplotlib figure directly into the tkinter window
        cv = FigureCanvasTkAgg(fig, master=p)
        cv.draw(); cv.get_tk_widget().pack(fill="both", expand=True, pady=4)

    # MODE 1 — Predictor by Mo
    def show_predictor(self):
        self.switch_to("Predictor")
        f = self.scrollable(self.content)

        tk.Label(f, text="Predictor  generate spectra for a known molecule",
                 bg=BG, fg=ACCENT, font=("Arial", 11, "bold")).pack(anchor="w", pady=(10,2))
        tk.Label(f, text="Note: only works for the 15 built-in database molecules.",
                 bg=BG, fg=GOLD, font=("Arial", 8, "italic")).pack(anchor="w", pady=(0,6))

        tk.Label(f, text="Molecule name:", bg=BG, fg=TEXT, font=("Arial", 10)).pack(anchor="w")
        name_entry = tk.Entry(f, width=30, bg=CARD, fg=WHITE, insertbackground=WHITE,
                              font=("Arial", 10), relief="flat")
        name_entry.pack(anchor="w", pady=2)

        tk.Label(f, text="Spectra to generate:", bg=BG, fg=TEXT,
                 font=("Arial", 10)).pack(anchor="w", pady=(8,2))
        choice = tk.StringVar(value="all")
        for v, l in [("all","All four"),("ir","IR"),("ms","MS"),("uv","UV-Vis"),("nmr","NMR")]:
            tk.Radiobutton(f, text=l, variable=choice, value=v, bg=BG, fg=TEXT,
                           selectcolor=CARD, activebackground=BG,
                           font=("Arial", 9)).pack(anchor="w")

        result_frame = tk.Frame(f, bg=BG)
        result_frame.pack(fill="both", expand=True, pady=6)

        def generate():
            name = name_entry.get().strip().lower()
            if not name:
                messagebox.showerror("Error", "Enter a molecule name.")
                return
            mol = self.db.get(name)
            if mol is None:
                messagebox.showerror("Not Found",
                    f"'{name}' not in database.\nAvailable: {', '.join(self.db.list_names())}")
                return
            for widget in result_frame.winfo_children():
                widget.destroy()
            plt.close("all")

            tk.Label(result_frame, text=mol.summary(), bg=BG, fg=GREEN,
                     font=("Courier", 9), justify="left").pack(anchor="w", pady=4)

            struct_path = f"output/structures/{mol.name.replace(' ','_')}.png"
            draw_molecule_2d(mol, save_path=struct_path, size=(300, 250))
            if os.path.exists(struct_path):
                img = Image.open(struct_path).resize((300, 250))
                tk_img = ImageTk.PhotoImage(img)
                img_label = tk.Label(result_frame, image=tk_img, bg=CARD)
                img_label.image = tk_img  # keep reference so image is not garbage collected
                img_label.pack(anchor="w", pady=6)

            loading = tk.Label(result_frame, text="Generating spectra...", bg=BG,
                               fg=ACCENT, font=("Arial", 10, "italic"))
            loading.pack(anchor="w", pady=4)

            selected = choice.get()

            # generate and embed spectra 
            loading.destroy()
            for key, Cls in [("ir", IRSpectrum), ("ms", MassSpectrum),
                              ("uv", UVVisSpectrum), ("nmr", NMRSpectrum)]:
                if selected in ("all", key):
                    fig = Cls(mol).plot()
                    if fig:
                        fig.savefig(f"output/spectra/{mol.name.replace(' ','_')}_{key.upper()}.png",
                                    dpi=120, bbox_inches="tight")
                        canvas = FigureCanvasTkAgg(fig, master=result_frame)
                        canvas.draw()
                        canvas.get_tk_widget().pack(fill="both", expand=True, pady=4)

        self.make_btn(f, "Generate", generate)

    """MODE 2 — Viewer (code written by Yoyo)"""
    def show_viewer(self):
        self.switch_to("Viewer")
        f = self.scrollable(self.content)

        tk.Label(f, text="Viewer draw and export any molecule structure",
                 bg=BG, fg=ACCENT, font=("Arial", 11, "bold")).pack(anchor="w", pady=(10,2))
        tk.Label(f, text="Enter any molecule name (database or online):",
                 bg=BG, fg=TEXT, font=("Arial", 10)).pack(anchor="w")
        name_entry = tk.Entry(f, width=40, bg=CARD, fg=WHITE, insertbackground=WHITE,
                              font=("Arial", 10), relief="flat")
        name_entry.pack(anchor="w", pady=2)

        info_label   = tk.Label(f, text="", bg=BG, fg=GREEN, font=("Arial", 9),
                                wraplength=520, justify="left")
        info_label.pack(anchor="w", pady=2)
        image_label  = tk.Label(f, bg=BG)
        image_label.pack(anchor="w", pady=4)
        status_label = tk.Label(f, text="", bg=BG, fg="#777777", font=("Arial", 8))
        status_label.pack(anchor="w")

        def draw():
            name = name_entry.get().strip().lower()
            if not name:
                messagebox.showerror("Error", "Enter a molecule name.")
                return
            mol, full = self.db.get_or_fetch(name)
            if mol is None:
                messagebox.showerror("Not Found",
                    f"Could not find '{name}'.\nCheck spelling or internet connection.")
                return
            self.viewer_mol  = mol
            self.viewer_name = name
            path = f"output/structures/{name.replace(' ','_')}.png"
            draw_molecule_2d(mol, save_path=path, size=(360, 300))
            if os.path.exists(path):
                img = Image.open(path).resize((360, 300))
                tk_img = ImageTk.PhotoImage(img)
                image_label.config(image=tk_img, bg=CARD)
                image_label.image = tk_img
            if full:
                info_label.config(text=(
                    f"{mol.iupac_name}   {mol.formula}   MW={mol.molecular_weight} g/mol   {mol.compound_class}\n"
                    f"Groups: {', '.join(mol.functional_groups)}\n"
                    f"IR peaks: {mol.ir_peaks}\n"
                    f"MS peaks: {mol.ms_peaks}\n"
                    f"UV peaks: {mol.uv_peaks}\n"
                    f"NMR peaks: {mol.nmr_peaks}\n"))
            else:
                info_label.config(text="Molecule not in database. Structure retrieved online. No spectral data available.")
            status_label.config(text=f"Image saved to {path}")

        """
        choice for user to export file (.xyz/.mol) to look the 3D molecule in Avogadro
        the user can choose from both
        """
        def export(filetype):
            if self.viewer_mol is None:
                messagebox.showerror("Error", "Draw a molecule first.")
                return
            from rdkit import Chem
            from rdkit.Chem import AllChem
            m3d = Chem.AddHs(self.viewer_mol.get_rdkit_mol())
            if AllChem.EmbedMolecule(m3d, AllChem.ETKDGv3()) == -1:
                messagebox.showerror("Error", "Could not generate 3D coordinates.")
                return
            AllChem.UFFOptimizeMolecule(m3d)
            fn = f"output/structures/{self.viewer_name.replace(' ','_')}.{filetype}"
            if filetype == "mol":
                Chem.MolToMolFile(m3d, fn)
            else:
                Chem.MolToXYZFile(m3d, fn)
            messagebox.showinfo("Exported", f"Saved as {fn}")
            status_label.config(text=f".{filetype} saved to {fn}")

        name_entry.bind("<Return>", lambda e: draw())
        btn_row = tk.Frame(f, bg=BG)
        btn_row.pack(anchor="w", pady=4)
        
        """button for export and draw"""
        self.make_btn(btn_row, "Draw", draw, side="left")
        self.make_btn(btn_row, "Export .mol for Avogadro", lambda: export("mol"), side="left")
        self.make_btn(btn_row, "Export .xyz for Avogadro", lambda: export("xyz"), side="left")
"""end of code written by Yoyo"""

    # MODE 3 — Database by Mo (just shows what's in the database tbh and a fun fact about it)
    def show_database(self):
        self.switch_to("Database")
        f = self.scrollable(self.content)

        tk.Label(f, text="Database  15 built-in compounds", bg=BG, fg=ACCENT,
                 font=("Arial", 11, "bold")).pack(anchor="w", pady=(10,6))

        for mol in sorted(self.db.all_molecules(), key=lambda m: m.name):
            card = tk.Frame(f, bg=CARD, pady=6, padx=12)
            card.pack(fill="x", pady=3)
            tk.Label(card, text=f"{mol.name.title()}  {mol.formula}",
                     bg=CARD, fg=ACCENT, font=("Arial", 10, "bold")).pack(anchor="w")
            if mol.fun_fact:
                tk.Label(card, text=f"Fun Fact: {mol.fun_fact}", bg=CARD, fg=GOLD,
                         font=("Arial", 8, "italic"), wraplength=700,
                         justify="left").pack(anchor="w")
                
""" MODE 4 — Comparator (written by Yoyo)"""
    def show_comparator(self):
        self.switch_to("Comparator")
        f = self.scrollable(self.content)
        tk.Label(f, text="Comparator compare two molecules side by side",
                 bg=BG, fg=ACCENT, font=("Arial", 11, "bold")).pack(anchor="w", pady=(10,2))
        tk.Label(f, text="Only able to compare the 15 molecules from database", bg=BG, fg=GOLD, font=("Arial", 8, "italic")).pack(anchor="w", pady=(0,4))

        row = tk.Frame(f, bg=BG); row.pack(anchor="w", pady=4)
        tk.Label(row, text="Molecule A:", bg=BG, fg=TEXT, font=("Arial", 10)).pack(side="left", padx=(0,0))
        a_entry = tk.Entry(row, width=22, bg=CARD, fg=WHITE, insertbackground=WHITE,
                           font=("Arial", 10), relief="flat"); a_entry.pack(side="left", padx=4)
        tk.Label(row, text="Molecule B:", bg=BG, fg=TEXT, font=("Arial", 10)).pack(side="left", padx=(10,0))
        b_entry = tk.Entry(row, width=22, bg=CARD, fg=WHITE, insertbackground=WHITE,
                           font=("Arial", 10), relief="flat"); b_entry.pack(side="left", padx=4)
        rf = tk.Frame(f, bg=BG); rf.pack(fill="both", expand=True, pady=6)

        def run():
            na, nb = a_entry.get().strip().lower(), b_entry.get().strip().lower() 
            """validate the inputs"""
            if not na or not nb: return messagebox.showerror("Error", "Enter both molecule names.") 
            """close previous plot"""
            for w in rf.winfo_children(): 
                w.destroy()
            plt.close("all")
            """define heavy work"""
            def heavy(): 
                return self.comparator.compare(na, nb)
            """show comparison"""
            def done(res):  
                fig, summary = res
                if fig is None: return messagebox.showerror("Error", summary)
                self.embed_fig(rf, fig)
                fig.savefig(f"output/spectra/comparison_{na}_vs_{nb}.png", dpi=120, bbox_inches="tight")
                tk.Label(rf, text="Summary", bg=BG, fg=ACCENT,
                         font=("Arial", 10, "bold")).pack(anchor="w", pady=(8,2))
                box = tk.Text(rf, bg=CARD, fg=TEXT, font=("Courier", 9), relief="flat",
                              height=4, width=62, wrap="word",
                              highlightthickness=1, highlightbackground=ACCENT)
                box.insert("1.0", summary); box.config(state="disabled")
                box.pack(anchor="w", pady=4)
            result = heavy()
            done(result)
        self.make_btn(f, "Compare", run)


if __name__ == "__main__":
    root = tk.Tk()
    SpectroscopyApp(root)
    root.mainloop()
