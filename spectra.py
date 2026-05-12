# Code by Mo
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from molecules import Molecule

# colour palette for graphs
IR_COLOR    = "#00BFFF"
MS_COLOR    = "#FF6B6B"
UV_COLOR    = "#9B59B6"
NMR_COLOR   = "#2ECC71"
AXIS_COLOR  = "#CCCCCC"
BG_COLOR    = "#1A1A2E"
PANEL_COLOR = "#16213E"

# dark background 
def apply_dark_style(ax, fig):
    # consistent dark theme across all spectrum plots
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor(PANEL_COLOR)
    ax.tick_params(colors=AXIS_COLOR, labelsize=9)
    ax.xaxis.label.set_color(AXIS_COLOR)
    ax.yaxis.label.set_color(AXIS_COLOR)
    ax.title.set_color("white")
    for spine in ax.spines.values():
        spine.set_edgecolor(AXIS_COLOR)
    ax.grid(True, linestyle="--", alpha=0.25, color=AXIS_COLOR)

# plotting the IR graph using Lorentzian broadening (will be referenced)
class IRSpectrum:
    
    FWHM = {"O-H": 200, "N-H": 80, "C-H": 25, "=C-H": 25,
            "C=O": 35, "C=C": 30, "C-O": 40, "default": 30}

    def __init__(self, molecule: Molecule):
        self.molecule = molecule

    def _get_fwhm(self, assignment):
        for key, fwhm in self.FWHM.items():
            if key in assignment:
                return float(fwhm)
        return float(self.FWHM["default"])

    def _build_spectrum(self, x):
        # sum Lorentzian contributions from each peak
        y = np.zeros_like(x)
        for wn, intensity, assignment in self.molecule.ir_peaks:
            g = self._get_fwhm(assignment)
            y += intensity * (g/2)**2 / ((x - wn)**2 + (g/2)**2)
        if y.max() > 0:
            y = y / y.max() * 100
        return y
#returning the plot
    def plot(self, fig=None, ax=None):
        # returns a matplotlib figure — caller embeds it in the GUI
        if not self.molecule.ir_peaks:
            return None
        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=(7, 3.5))
        apply_dark_style(ax, fig)
        x = np.linspace(4000, 500, 4000)
        # IR shown as transmittance — convention is peaks go down
        transmittance = 100 - self._build_spectrum(x)
        ax.plot(x, transmittance, color=IR_COLOR, linewidth=1.5)
        for wn, intensity, _ in self.molecule.ir_peaks:
            t = 100 - intensity
            ax.axvline(x=wn, ymin=0, ymax=t/100, color=IR_COLOR, alpha=0.4, linestyle="--", linewidth=0.8)
            ax.annotate(f"{wn:.0f}", xy=(wn, t), xytext=(0, -18),
                        textcoords="offset points", ha="center", fontsize=6, color=IR_COLOR, rotation=90)
        # x-axis right to left — standard IR convention
        ax.set_xlim(4000, 500)
        ax.invert_xaxis()
        ax.set_ylim(-5, 110)
        ax.set_xlabel("Wavenumber (cm-1)", fontsize=9)
        ax.set_ylabel("% Transmittance", fontsize=9)
        ax.set_title(f"IR — {self.molecule.name.title()}", fontsize=10, fontweight="bold")
        ax.axvspan(3600, 2500, alpha=0.05, color="yellow", label="O-H/N-H/C-H")
        ax.axvspan(2000, 1500, alpha=0.05, color="cyan",   label="Carbonyl/C=C")
        ax.axvspan(1500,  400, alpha=0.05, color="green",  label="Fingerprint")
        ax.legend(loc="upper left", fontsize=7, facecolor=PANEL_COLOR, labelcolor=AXIS_COLOR)
        plt.tight_layout()
        return fig

# mass spec plot
class MassSpectrum:
    def __init__(self, molecule: Molecule):
        self.molecule = molecule

    def plot(self, fig=None, ax=None):
        if not self.molecule.ms_peaks:
            return None
        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=(7, 3.5))
        apply_dark_style(ax, fig)
        mz_values   = [p[0] for p in self.molecule.ms_peaks]
        intensities = [p[1] for p in self.molecule.ms_peaks]
        base_peak   = self.molecule.get_base_peak()
        ax.vlines(mz_values, 0, intensities, color=MS_COLOR, linewidth=3.5, alpha=0.85)
        ax.scatter(mz_values, intensities, color=MS_COLOR, zorder=5, s=20)
        for mz, intensity, label in zip(mz_values, intensities, [p[2] for p in self.molecule.ms_peaks]):
            # highlight base peak in gold
            colour = "#FFD700" if (base_peak and abs(mz - base_peak[0]) < 0.5) else AXIS_COLOR
            ax.annotate(f"m/z {mz:.0f}", xy=(mz, intensity), xytext=(0, 6),
                        textcoords="offset points", ha="center", fontsize=7, color=colour)
        mol_ion = self.molecule.get_molecular_ion()
        if mol_ion:
            mol_i = next((p[1] for p in self.molecule.ms_peaks if abs(p[0] - mol_ion) < 0.5), 0)
            ax.annotate("M+", xy=(mol_ion, mol_i), xytext=(mol_ion+3, mol_i+8),
                        fontsize=9, color="#FFD700",
                        arrowprops=dict(arrowstyle="->", color="#FFD700", lw=1.2))
        ax.set_xlim(0, max(mz_values)*1.15)
        ax.set_ylim(0, 115)
        ax.set_xlabel("m/z", fontsize=9)
        ax.set_ylabel("Relative Abundance (%)", fontsize=9)
        ax.set_title(f"MS — {self.molecule.name.title()}", fontsize=10, fontweight="bold")
        plt.tight_layout()
        return fig

#uv-vis plot 
class UVVisSpectrum:
    # Gaussian broadening — correct for electronic transitions in solution, see notes.md
    SIGMA = {"pi->pi*": 15, "n->pi*": 10, "n->sigma*": 12, "default": 20}

    def __init__(self, molecule: Molecule):
        self.molecule = molecule

    def _get_sigma(self, transition):
        for key, s in self.SIGMA.items():
            if key in transition:
                return float(s)
        return float(self.SIGMA["default"])

    def _build_spectrum(self, wl):
        abs_ = np.zeros_like(wl)
        for lam0, eps, transition in self.molecule.uv_peaks:
            s = self._get_sigma(transition)
            # log10 epsilon keeps scale manageable across weak and strong bands
            abs_ += np.log10(max(eps, 1)) * np.exp(-(wl - lam0)**2 / (2*s**2))
        return abs_

    def plot(self, fig=None, ax=None):
        if not self.molecule.uv_peaks:
            return None
        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=(7, 3.5))
        apply_dark_style(ax, fig)
        wl   = np.linspace(200, 800, 1200)
        abs_ = self._build_spectrum(wl)
        ax.plot(wl, abs_, color=UV_COLOR, linewidth=2.0)
        ax.fill_between(wl, 0, abs_, color=UV_COLOR, alpha=0.15)
        for lam0, eps, transition in self.molecule.uv_peaks:
            a = np.log10(max(eps, 1))
            ax.annotate(f"l={lam0:.0f}nm\n({transition})", xy=(lam0, a),
                        xytext=(lam0+15, a+0.2), fontsize=7, color=UV_COLOR,
                        arrowprops=dict(arrowstyle="->", color=UV_COLOR, lw=0.8))
        # shade visible light region
        for lo, hi, c in [(380,450,"violet"),(450,495,"blue"),(495,570,"green"),(570,620,"yellow"),(620,750,"red")]:
            ax.axvspan(lo, hi, alpha=0.07, color=c)
        ax.set_xlim(200, 800)
        ax.set_ylim(0, max(abs_)*1.4+0.1)
        ax.set_xlabel("Wavelength (nm)", fontsize=9)
        ax.set_ylabel("Absorbance (log10 e)", fontsize=9)
        ax.set_title(f"UV-Vis — {self.molecule.name.title()}", fontsize=10, fontweight="bold")
        plt.tight_layout()
        return fig

# NMR plot
class NMRSpectrum:
    # first-order multiplet splitting using Pascal's triangle, J=7Hz at 300MHz, see notes.md
    SPECTROMETER_MHZ = 300
    J_PPM  = 7.0 / 300
    FWHM_NMR = 0.008

    def __init__(self, molecule: Molecule):
        self.molecule = molecule
        J = self.J_PPM
        # offsets and Pascal coefficients per multiplicity
        self._mult = {
            "s":    ([0.0],                          [1]),
            "d":    ([-J/2, J/2],                    [1,1]),
            "t":    ([-J, 0, J],                      [1,2,1]),
            "q":    ([-1.5*J,-0.5*J,0.5*J,1.5*J],   [1,3,3,1]),
            "dd":   ([-J, 0, J],                      [1,2,1]),
            "dq":   ([-J, 0, J],                      [1,2,1]),
            "ddt":  ([-J, 0, J],                      [1,2,1]),
            "sept": ([-3*J,-2*J,-J,0,J,2*J,3*J],     [1,6,15,20,15,6,1]),
            "m":    ([0.0],                           [1]),
        }

    def _lorentzian(self, x, centre, amplitude, fwhm):
        return amplitude * (fwhm/2)**2 / ((x - centre)**2 + (fwhm/2)**2)

    def _build_spectrum(self, ppm_axis):
        y = np.zeros_like(ppm_axis)
        for shift, mult, integration, _ in self.molecule.nmr_peaks:
            offsets, coeffs = self._mult.get(mult.lower(), ([0.0],[1]))
            fwhm = self.FWHM_NMR * (3 if mult.lower() == "m" else 1)
            coeff_sum = sum(coeffs)
            for offset, coeff in zip(offsets, coeffs):
                y += self._lorentzian(ppm_axis, shift+offset, integration*(coeff/coeff_sum), fwhm)
        return y

    def plot(self, fig=None, ax=None):
        if not self.molecule.nmr_peaks:
            return None
        if fig is None or ax is None:
            fig, ax = plt.subplots(figsize=(7, 3.5))
        apply_dark_style(ax, fig)
        shifts  = [p[0] for p in self.molecule.nmr_peaks]
        ppm_max = max(shifts) + 0.8
        ppm_min = max(min(shifts) - 0.5, 0)
        ppm_axis = np.linspace(ppm_min, ppm_max+1.0, 5000)
        y = self._build_spectrum(ppm_axis)
        ax.plot(ppm_axis, y, color=NMR_COLOR, linewidth=1.8)
        ax.fill_between(ppm_axis, 0, y, alpha=0.2, color=NMR_COLOR)
        for shift, mult, integ, assign in self.molecule.nmr_peaks:
            peak_y = self._build_spectrum(np.array([shift]))[0]
            ax.annotate(f"d{shift:.2f}\n({mult},{integ:.0f}H)\n{assign}",
                        xy=(shift, peak_y), xytext=(shift, peak_y+0.3), ha="center",
                        fontsize=6, color=NMR_COLOR,
                        bbox=dict(boxstyle="round,pad=0.2", facecolor=PANEL_COLOR, alpha=0.7))
        # NMR x-axis runs right to left, TMS at 0
        ax.set_xlim(ppm_max+0.5, ppm_min-0.2)
        ax.invert_xaxis()
        ax.set_xlabel("Chemical Shift d (ppm)", fontsize=9)
        ax.set_ylabel("Intensity", fontsize=9)
        ax.set_title(f"1H-NMR — {self.molecule.name.title()} (300 MHz)", fontsize=10, fontweight="bold")
        plt.tight_layout()
        return fig
