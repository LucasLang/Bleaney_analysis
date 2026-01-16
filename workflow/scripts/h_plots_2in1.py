import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
from py_functions import plot_error_vs_h

h_minexp = int(sys.argv[1])
h_maxexp = int(sys.argv[2])
h_norms_file_TM = sys.argv[3]
h_norms_file_Ln = sys.argv[4]
h_plot_file = sys.argv[5]

# Define the results directory
results_dir = "results"

h_vals = [10**h for h in range(h_minexp, h_maxexp + 1)]
# Plotting
fig, ax = plt.subplots(figsize=(5,3))
fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(6.7, 3)
)

plot_error_vs_h(axes[0], h_vals, h_norms_file_TM, legend=False)    # it is sufficient to have the legend in one of the two plots
plot_error_vs_h(axes[1], h_vals, h_norms_file_Ln)

label_size = mpl.rcParams['axes.labelsize']
axes[0].set_title(r"$S=0$, $D=3$ cm$^{-1}$, $E/D=0$", fontsize=label_size)
axes[1].set_title("Yb(III) complex", fontsize=label_size)


plt.tight_layout()

# Save the plot as a PNG file in the results directory
plt.savefig(h_plot_file, dpi=300, bbox_inches='tight')
