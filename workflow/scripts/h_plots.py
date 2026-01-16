import matplotlib.pyplot as plt
import numpy as np
import sys
from py_functions import plot_error_vs_h

h_minexp = int(sys.argv[1])
h_maxexp = int(sys.argv[2])
h_norms_file = sys.argv[3]
h_plot_file = sys.argv[4]

# Define the results directory
results_dir = "results"

h_vals = [10**h for h in range(h_minexp, h_maxexp + 1)]
# Plotting
fig, ax = plt.subplots(figsize=(5,3))

plot_error_vs_h(ax, h_vals, h_norms_file)


plt.tight_layout()

# Save the plot as a PNG file in the results directory
plt.savefig(h_plot_file, dpi=300, bbox_inches='tight')
