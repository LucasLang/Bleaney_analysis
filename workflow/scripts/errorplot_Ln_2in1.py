import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.offsetbox import AnchoredText
import sys
from py_functions import convergence_plot

Tmin = int(sys.argv[1])
Tmax = int(sys.argv[2])
Tinterval = int(sys.argv[3])
diff_norms_files = [sys.argv[4], sys.argv[6]]
lns = [sys.argv[5], sys.argv[7]]
error_analysis_plot_file = sys.argv[8]

temps = np.arange(Tmin, Tmax+1, Tinterval)

# --- load data (6 columns: beta^1..beta^6) ---
data_matrices = [np.loadtxt(file) for file in diff_norms_files]

# --- axis labels for x ---
beta_terms = [r'$\beta$', r'$\beta^2$', r'$\beta^3$', r'$\beta^4$', r'$\beta^5$', r'$\beta^6$']

# --- colors by temperature (log-normalized) ---
norm = mcolors.Normalize(vmin=min(temps), vmax=max(temps))

fig = plt.figure(figsize=(6.7, 3), layout="constrained")
gs = fig.add_gridspec(nrows=2, ncols=2, height_ratios=[0.18, 1.0])

axes = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
ax_legend = fig.add_subplot(gs[0, :])
ax_legend.axis("off")

formatted_labels = [f'{T:.0f}' for T in temps]
for i in [0, 1]:
    convergence_plot(axes[i], beta_terms, data_matrices[i], formatted_labels, temps, norm)

    axes[i].set_ylim(1e-5, 1e1)
    plotlabels = f"{lns[i]}(III)"
    at = AnchoredText(
        plotlabels,
        loc="lower left",
        frameon=True,
        borderpad=0.6,
        pad=0.4,
    )

    axes[i].add_artist(at)

# --- temperature legend (two columns) ---
handles, labels = axes[0].get_legend_handles_labels()     # same labels for each Ln ion
ax_legend.legend(
    handles, labels,
    title=r"$T$/K",
    loc='center',
    frameon=True,
    ncol = 5
)

plt.tight_layout()

# --- save ---
plt.savefig(error_analysis_plot_file, dpi=300, bbox_inches='tight')

