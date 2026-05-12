# By Yoyo
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from molecules import Molecule, MoleculeDatabase
from spectra import (IRSpectrum, MassSpectrum, UVVisSpectrum, NMRSpectrum,
                     apply_dark_style, BG_COLOR, PANEL_COLOR, AXIS_COLOR,
                     IR_COLOR, MS_COLOR, UV_COLOR, NMR_COLOR)

# tolerance for IR peak comparison between two molecules
IR_COMPARE_TOL = 50


class SpectrumComparator:
    def __init__(self, database: MoleculeDatabase):
        self.db = database

    def compare(self, name_a, name_b):
        """
        returns (figure, summary_text) or error message
        """
        mol_a = self.db.get(name_a)
        mol_b = self.db.get(name_b)
        if mol_a is None: return None, f"'{name_a}' not found in database."
        if mol_b is None: return None, f"'{name_b}' not found in database."
        fig = self._build_figure(mol_a, mol_b)
        return fig, self._build_summary(mol_a, mol_b)

    def _build_figure(self, mol_a, mol_b):
      """
      shows spectra in an organised layout
      """
        fig = plt.figure(figsize=(14, 16))
        fig.patch.set_facecolor(BG_COLOR)
        # GridSpec for subplot layout
        gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.55, wspace=0.3, top=0.93, bottom=0.04)

        # store plot functions in list to avoid repeating the loop 4 times
        # heavy work with loop may get crash
        
        plotters = [
            (IR_COLOR, "IR", self._plot_ir),
            (MS_COLOR, "MS", self._plot_ms),
            (UV_COLOR, "UV-Vis", self._plot_uv),
            (NMR_COLOR, "NMR", self._plot_nmr),
        ]
        for row, (colour, label, fn) in enumerate(plotters):
            ax_a = fig.add_subplot(gs[row, 0])
            ax_b = fig.add_subplot(gs[row, 1])
            apply_dark_style(ax_a, fig)
            apply_dark_style(ax_b, fig)
            fn(ax_a, mol_a, colour, "A")
            fn(ax_b, mol_b, colour, "B")
            ax_a.set_ylabel(label, fontsize=9, color=colour, fontweight="bold")

        fig.suptitle(
            f"A: {mol_a.name.title()} ({mol_a.formula})   vs   B: {mol_b.name.title()} ({mol_b.formula})",
            fontsize=13, fontweight="bold", color="white", y=0.97
        )
        return fig

    def _no_data(self, ax, mol, colour, label, msg):
        ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes, color=AXIS_COLOR)
        ax.set_title(f"{label}: {mol.name.title()}", fontsize=9, color=colour)

    
    def _plot_ir(self, ax, mol, colour, label):
      """
      plot IR spectrum
      """
        if not mol.ir_peaks: return self._no_data(ax, mol, colour, label, "No IR data")
        x = np.linspace(4000, 500, 4000)
        ax.plot(x, 100 - IRSpectrum(mol)._build_spectrum(x), color=colour, linewidth=1.3)
        ax.set_xlim(4000, 500); ax.invert_xaxis(); ax.set_ylim(-5, 110)
        ax.set_xlabel("Wavenumber (cm-1)", fontsize=7)
        ax.set_title(f"{label}: {mol.name.title()}", fontsize=8, color=colour)
        for wn, inten, _ in mol.ir_peaks:
            ax.plot(wn, 100-inten, "o", color=colour, markersize=2)

    def _plot_ms(self, ax, mol, colour, label):
      """
      plot MS spectrum
      """
        if not mol.ms_peaks: return self._no_data(ax, mol, colour, label, "No MS data")
        mz  = [p[0] for p in mol.ms_peaks]
        rel = [p[1] for p in mol.ms_peaks]
        ax.vlines(mz, 0, rel, color=colour, linewidth=2.5, alpha=0.85)
        ax.scatter(mz, rel, color=colour, s=10, zorder=5)
        for m, r in zip(mz, rel):
            ax.annotate(f"{m:.0f}", xy=(m, r), xytext=(0,4), textcoords="offset points",
                        fontsize=6, ha="center", color=colour)
        ax.set_xlim(0, max(mz)*1.15); ax.set_ylim(0, 115)
        ax.set_xlabel("m/z", fontsize=7)
        ax.set_title(f"{label}: {mol.name.title()}", fontsize=8, color=colour)

    def _plot_uv(self, ax, mol, colour, label):
      """
      plot UV spectrum
      """
        if not mol.uv_peaks: return self._no_data(ax, mol, colour, label, "No UV data")
        wl   = np.linspace(200, 800, 1200)
        abs_ = UVVisSpectrum(mol)._build_spectrum(wl)
        ax.plot(wl, abs_, color=colour, linewidth=1.3)
        ax.fill_between(wl, 0, abs_, color=colour, alpha=0.12)
        ax.set_xlim(200, 800); ax.set_xlabel("Wavelength (nm)", fontsize=7)
        ax.set_title(f"{label}: {mol.name.title()}", fontsize=8, color=colour)

    def _plot_nmr(self, ax, mol, colour, label):
      """
      plot NMR spectrum
      """
        if not mol.nmr_peaks: return self._no_data(ax, mol, colour, label, "No NMR data")
        shifts  = [p[0] for p in mol.nmr_peaks]
        ppm_max = max(shifts)+0.8; ppm_min = max(min(shifts)-0.5, 0)
        ppm_axis = np.linspace(ppm_min, ppm_max+1.0, 3000)
        y = NMRSpectrum(mol)._build_spectrum(ppm_axis)
        ax.plot(ppm_axis, y, color=colour, linewidth=1.3)
        ax.fill_between(ppm_axis, 0, y, color=colour, alpha=0.12)
        ax.set_xlim(ppm_max+0.5, ppm_min-0.2); ax.invert_xaxis()
        ax.set_xlabel("d (ppm)", fontsize=7)
        ax.set_title(f"{label}: {mol.name.title()}", fontsize=8, color=colour)

    def _build_summary(self, mol_a, mol_b):
        """
        generates the difference summary in text and show in the GUI
        """
        lines = []

        if mol_a.formula == mol_b.formula:
            lines.append(f"Same formula ({mol_a.formula}) — structural isomers.")
        else:
            lines.append(f"Formula: A={mol_a.formula}  B={mol_b.formula}")
            lines.append(f"MW: A={mol_a.molecular_weight}  B={mol_b.molecular_weight} g/mol")

        only_a = set(mol_a.functional_groups) - set(mol_b.functional_groups)
        only_b = set(mol_b.functional_groups) - set(mol_a.functional_groups)
        common = set(mol_a.functional_groups) & set(mol_b.functional_groups)
        if common:  lines.append(f"Shared groups: {', '.join(sorted(common))}")
        if only_a:  lines.append(f"Only in A: {', '.join(sorted(only_a))}")
        if only_b:  lines.append(f"Only in B: {', '.join(sorted(only_b))}")

        
        return "\n".join(lines)
