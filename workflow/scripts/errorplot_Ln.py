import matplotlib.pyplot as plt
import numpy as np
import json
import os
import re
import glob
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# --- detect lanthanoid name from file like data/Bkq_<LN>_real ---
lnfile = glob.glob("data/Bkq_*_real")[0]
lnname = os.path.basename(lnfile)
match = re.search(r"Bkq_(\w+)_real", lnname)
ln = match.group(1) if match else "Unknown"

# --- results dir and temps from runs.json ---
results_dir = "results"
runs_file = os.path.join(results_dir, "runs.json")
with open(runs_file, "r") as file:
    runs = json.load(file)

temps = [T[0] for T in runs]  # list of temperatures (K)

# --- load data (6 columns: beta^1..beta^6) ---
file_paths = [os.path.join(results_dir, f"diff_norm{i}.txt") for i in range(1, 7)]
data_matrix = [np.loadtxt(path) for path in file_paths]
data_matrix = np.array(data_matrix).T  # shape (n_curves, 6)

# --- axis labels for x ---
beta_terms = [r'$\beta$', r'$\beta^2$', r'$\beta^3$', r'$\beta^4$', r'$\beta^5$', r'$\beta^6$']

# --- colors by temperature (log-normalized) ---
norm = mcolors.LogNorm(vmin=min(temps), vmax=max(temps))
cmap = cm.viridis
colors = [cmap(norm(t)) for t in temps]

# --- figure/axes ---
fig, ax = plt.subplots(figsize=(10, 6))

# --- plot each temperature curve ---
for i, (temp_data, color) in enumerate(zip(data_matrix, colors)):
    ax.plot(beta_terms, temp_data, marker='o', color=color, label=f'{temps[i]:.0f}')

# --- axes formatting ---
ax.set_yscale('log')
ax.set_ylim(1e-5, 1e1)
ax.set_xlabel("Highest-order Term in the Expansion", fontsize=12, fontweight='bold')
ax.set_ylabel("Relative Error Values", fontsize=12, fontweight='bold')
ax.grid(True, which="major", ls="-", alpha=0.2)
ax.tick_params(axis='y', which='minor', length=0)

# --- temperature legend (two columns) ---
handles, labels = ax.get_legend_handles_labels()
temp_leg = ax.legend(
    handles, labels,
    title="Temperatures in Kelvin",
    ncol=2,
    loc='lower left',
    bbox_to_anchor=(0.02, 0.08),
    fontsize=10,
    title_fontsize=11,
    frameon=True,
    fancybox=True,
    shadow=True,
    columnspacing=1.2,
    handlelength=1.8,
    handletextpad=0.6,
    borderaxespad=0.0,
    labelspacing=0.5
)
ax.add_artist(temp_leg)

# --- tidy header ABOVE the legend (robust text, not a second legend) ---
# Align left edges by using the same x in axes coordinates; lift y a bit above the legend
ax.text(
    0.05, 0.35, rf"Lanthanoid: {ln}",
    transform=ax.transAxes,
    fontsize=12,
    fontweight='bold',
    va='bottom', ha='left'
)

plt.tight_layout()

# --- save ---
error_analysis_plot_file = os.path.join(results_dir, "error_analysis_plot.png")
plt.savefig(error_analysis_plot_file, dpi=300, bbox_inches='tight')

