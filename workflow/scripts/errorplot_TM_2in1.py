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

S = [float(sys.argv[1]), float(sys.argv[6])]
E = [int(sys.argv[2]), int(sys.argv[7])]
D = [int(sys.argv[3]), int(sys.argv[8])]
diff_norms_files = [sys.argv[4], sys.argv[9]]
Tpointsfiles = [sys.argv[5], sys.argv[10]]

error_analysis_plot_file = sys.argv[11]

temps = [np.loadtxt(file) for file in Tpointsfiles]

data_matrices = [np.loadtxt(file) for file in diff_norms_files]

# X-axis labels
beta_terms = [r'$\beta$', r'$\beta^2$', r'$\beta^3$', r'$\beta^4$', r'$\beta^5$', r'$\beta^6$']

# Create figure and axes with larger size
fig = plt.figure(figsize=(6.7, 3), layout="constrained")
gs = fig.add_gridspec(nrows=1, ncols=3, width_ratios=[0.45, 0.45, 0.1])

axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
ax_legend = fig.add_subplot(gs[0, 2])
ax_legend.axis("off")

# --- colors by temperature (log-normalized) ---
for i in [0, 1]:
    norm = mcolors.LogNorm(vmin=min(temps[i]), vmax=max(temps[i]))
    temperature_ratios = [D[i] / (kB * T) for T in temps[i]]
    formatted_labels = [sci_label(float(T_ratio), sig=2) for T_ratio in temperature_ratios]
    convergence_plot(axes[i], beta_terms, data_matrices[i], formatted_labels, temps[i], norm)

    # Convert to fraction
    E_over_D_fraction = Fraction(E[i],D[i])

    twotimesS = int(2*S[i])
    S_frac = Fraction(twotimesS, 2)  # Fraction objects nicely interpolate into strings

    # Main plot formatting
    axes[i].set_ylim(1e-14, 1e6)
    plotlabels = fr"""$E/D = {E_over_D_fraction}$
$S = {S_frac}$"""

    at = AnchoredText(
        plotlabels,
        loc="lower left",
        frameon=True,
        borderpad=0.6,
        pad=0.4,
    )

    axes[i].add_artist(at)


# Global legend
handles, labels = axes[0].get_legend_handles_labels()     # We assume that axes have the same labels: do not mix positive and negative D!
ax_legend.legend(handles, labels, title=r"$D/k_\mathrm{B}T$",
          loc='center',
          #bbox_to_anchor=(0.01, 0.09),
          borderaxespad=0,
          frameon=True)


plt.tight_layout()

# Save with higher resolution and larger size
plt.savefig(error_analysis_plot_file, dpi=300, bbox_inches='tight')