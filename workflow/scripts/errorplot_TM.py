import matplotlib.pyplot as plt
import numpy as np
import json
import os
from fractions import Fraction

def sci_label(x, sig=2):
    """Return a MathText string like r'$m\times10^{e}$' for a number x."""
    if x == 0 or not np.isfinite(x):
        return r"$0$"
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / (10 ** exp)
    # Keep the sign with the mantissa
    return rf"${mant:.{sig}f}\times 10^{{{exp}}}$"

# Define the directories
results_dir = "results"
data_dir = "data"

ratio_file = os.path.join(data_dir, "ratio.txt")
with open(ratio_file, "r") as file:
    ratio = file.read()

k = 3.166811563e-6 
conversion_factor = 219474.63
kB = k * conversion_factor

# Read temperatures from runs.json in the results folder
runs_file = os.path.join(results_dir, "runs.json")
with open(runs_file, "r") as file:
    runs = json.load(file)

temps = [T[0] for T in runs]

D = np.loadtxt("data/D.txt", dtype=np.float64)
S = np.loadtxt("data/S.txt", dtype=np.float64)

# Define file paths in the 'results' folder
file_paths = [
    os.path.join(results_dir, "diff_norm1.txt"),
    os.path.join(results_dir, "diff_norm2.txt"),
    os.path.join(results_dir, "diff_norm3.txt"),
    os.path.join(results_dir, "diff_norm4.txt"),
    os.path.join(results_dir, "diff_norm5.txt"),
    os.path.join(results_dir, "diff_norm6.txt")
]

# Load data 
data_matrix = []
for path in file_paths:
    data = np.loadtxt(path) 
    data_matrix.append(data)

# Load data into matrices and transpose
data_matrix = np.array(data_matrix).T

# X-axis labels
beta_terms = [r'$\beta$', r'$\beta^2$', r'$\beta^3$', r'$\beta^4$', r'$\beta^5$', r'$\beta^6$']
temperature_points = temps

# Create figure and axes with larger size
plt.figure(figsize=(13, 10))

# Create two subplots: main plot and legend
#gs = plt.GridSpec(1, 2, width_ratios=[3, 1])
ax = plt.gca() #plt.subplot(gs[0])
#leg = plt.subplot(gs[1])
#leg.axis('off')

# Color map for better distinction between lines
colors = plt.cm.viridis(np.linspace(0, 1, len(data_matrix)))

# Plotting with improved visibility
for i, (temp_data, color) in enumerate(zip(data_matrix, colors)):
    temperature_ratio = D / (kB * temperature_points[i])
    formatted_label = sci_label(float(temperature_ratio), sig=2)
    ax.plot(beta_terms, temp_data, marker='o', label=formatted_label, 
            color=color, linewidth=2, markersize=8)

# Convert to fraction
E_over_D_fraction = Fraction(float(ratio)).limit_denominator()

S_frac = Fraction(float(S)).limit_denominator()

# Main plot formatting
ax.set_yscale('log')
ax.set_ylim(1e-14, 1e6)
ax.set_xlabel("Highest-order Term in the Expansion", fontsize=14, fontweight='bold')
ax.set_ylabel("Relative Error Values", fontsize=14, fontweight='bold')
#ax.set_title(f"TM", fontsize=14, fontweight='bold', pad=20)
ax.grid(True, which="both", ls="-", alpha=0.2)
ax.tick_params(axis='both', which='major', labelsize=15)

# First legend (D/kT values)
handles, labels = ax.get_legend_handles_labels()
first_legend = ax.legend(handles, labels, title="D/kT", 
          title_fontsize=14, fontsize=12,
          loc='lower left',
          bbox_to_anchor=(0.01, 0.09),
          borderaxespad=0,
          frameon=True,
          fancybox=True,
          shadow=True)



# Create dummy lines without actual lines for the second legend
dummy1 = plt.Line2D([], [], linestyle='none', label=f"E/D = {E_over_D_fraction}")
dummy2 = plt.Line2D([], [], linestyle='none', label=f"S = {S_frac}")

# Second legend (E/D as fraction and S)
secondary_legend = ax.legend(
    handles=[dummy1, dummy2],
    loc='lower left', 
    frameon=True, 
    fontsize=13
)

ax.add_artist(first_legend)  

plt.tight_layout()

# Save with higher resolution and larger size
error_analysis_plot_file = os.path.join(results_dir, "error_analysis_plot.png")
plt.savefig(error_analysis_plot_file, dpi=300, bbox_inches='tight')