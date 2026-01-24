import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.offsetbox import AnchoredText
import numpy as np
import sys
from fractions import Fraction
from py_functions import convergence_plot, sci_label

k = 3.166811563e-6 
conversion_factor = 219474.63
kB = k * conversion_factor

S = float(sys.argv[1])
E = int(sys.argv[2])
D = int(sys.argv[3])
diff_norms_file = sys.argv[4]
Tpointsfile = sys.argv[5]
error_analysis_plot_file = sys.argv[6]

temps = np.loadtxt(Tpointsfile)

data_matrix = np.loadtxt(diff_norms_file)

# X-axis labels
beta_terms = [r'$\beta$', r'$\beta^2$', r'$\beta^3$', r'$\beta^4$', r'$\beta^5$', r'$\beta^6$']

# Create figure and axes with larger size
fig = plt.figure(figsize=(5,3), layout = "constrained")
gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[0.8, 0.2])

ax = fig.add_subplot(gs[0, 0])
ax_legend = fig.add_subplot(gs[0, 1])
ax_legend.axis("off")

# --- colors by temperature (log-normalized) ---
norm = mcolors.LogNorm(vmin=min(temps), vmax=max(temps))
temperature_ratios = [D / (kB * T) for T in temps]
formatted_labels = [sci_label(float(T_ratio), sig=1) for T_ratio in temperature_ratios]
convergence_plot(ax, beta_terms, data_matrix, formatted_labels, temps, norm)

# Convert to fraction
E_over_D_fraction = Fraction(E,D)

twotimesS = int(2*S)
S_frac = Fraction(twotimesS, 2)  # Fraction objects nicely interpolate into strings

# Main plot formatting
ax.set_ylim(1e-14, 1e6)

# First legend (D/kT values)
handles, labels = ax.get_legend_handles_labels()
ax_legend.legend(handles, labels, title=r"$D/k_\mathrm{B}T$",
          loc='center',
          #bbox_to_anchor=(0.01, 0.09),
          borderaxespad=0,
          frameon=True)

plotlabels = fr"""$E/D = {E_over_D_fraction}$
$S = {S_frac}$"""

at = AnchoredText(
    plotlabels,
    loc="lower left",
    frameon=True,
    borderpad=0.6,
    pad=0.4,
)

ax.add_artist(at)

plt.tight_layout()

# Save with higher resolution and larger size
plt.savefig(error_analysis_plot_file, dpi=300, bbox_inches='tight')