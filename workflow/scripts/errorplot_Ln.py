import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.offsetbox import AnchoredText
import sys

Tmin = int(sys.argv[1])
Tmax = int(sys.argv[2])
Tinterval = int(sys.argv[3])
diff_norms_file = sys.argv[4]
ln = sys.argv[5]
error_analysis_plot_file = sys.argv[6]
temps = np.arange(Tmin, Tmax+1, Tinterval)

# --- load data (6 columns: beta^1..beta^6) ---
data_matrix = np.loadtxt(diff_norms_file)

# --- axis labels for x ---
beta_terms = [r'$\beta$', r'$\beta^2$', r'$\beta^3$', r'$\beta^4$', r'$\beta^5$', r'$\beta^6$']

# --- colors by temperature (log-normalized) ---
norm = mcolors.LogNorm(vmin=min(temps), vmax=max(temps))
cmap = cm.viridis
colors = [cmap(norm(t)) for t in temps]

fig = plt.figure(figsize=(5,3), layout = "constrained")
gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[0.8, 0.2])
ax = fig.add_subplot(gs[0, 0])
ax_legend = fig.add_subplot(gs[0, 1])
ax_legend.axis("off")

# --- plot each temperature curve ---
for i, (temp_data, color) in enumerate(zip(data_matrix, colors)):
    ax.plot(beta_terms, temp_data, marker='o', color=color, label=f'{temps[i]:.0f}')

# --- axes formatting ---
ax.set_yscale('log')
ax.set_ylim(1e-5, 1e1)
ax.set_xlabel("Highest-order term in the expansion")
ax.set_ylabel("Relative error")
ax.grid(True, which="major", ls="-", alpha=0.2)
ax.tick_params(axis='y', which='minor', length=0)

# --- temperature legend (two columns) ---
handles, labels = ax.get_legend_handles_labels()
ax_legend.legend(
    handles, labels,
    title=r"$T$/K",
    loc='center',
    frameon=True
)

plotlabels = f"{ln}(III)"
at = AnchoredText(
    plotlabels,
    loc="lower left",
    frameon=True,
    borderpad=0.6,
    pad=0.4,
)

ax.add_artist(at)

plt.tight_layout()

# --- save ---
plt.savefig(error_analysis_plot_file, dpi=300, bbox_inches='tight')

