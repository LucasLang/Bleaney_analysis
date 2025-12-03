import matplotlib.pyplot as plt
import numpy as np
import os
import sys

h_minexp = int(sys.argv[1])
h_maxexp = int(sys.argv[2])

# Define the results directory
results_dir = "results"

h_vals = [10**h for h in range(h_minexp, h_maxexp + 1)]

data_matrix = np.loadtxt(os.path.join(results_dir, "h_norms.txt"))   # Shape: hvals x 3

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))

beta_labels = [r'$\beta^1$', r'$\beta^2$', r'$\beta^3$'] 

# Plot each row in data_matrix (norms for each beta term)
Norders = data_matrix.shape[1]
for i in range(0,Norders):
    ax.plot(h_vals, data_matrix[:,i], marker='o', label=beta_labels[i])

# Set axes and labels
ax.set_xscale('log')  # Logarithmic scale for h values
ax.set_yscale('log')  # Logarithmic scale for norm values
ax.set_xlabel("h Values")
ax.set_ylabel("Relative Error Values")
#ax.set_title("Norm Values vs. h Terms")
ax.legend(title="Beta Terms", loc='best')

#plt.grid(visible=True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
#plt.show()

# Save the plot data to a text file in the results directory
plot_data_file = os.path.join(results_dir, "plot_data.txt")
np.savetxt(plot_data_file, np.column_stack((h_vals, data_matrix)), header='h_vals beta1 beta2 beta3')

# Save the plot as a PNG file in the results directory
h_plot_file = os.path.join(results_dir, "h_plot.png")
plt.savefig(h_plot_file, dpi=300, bbox_inches='tight')