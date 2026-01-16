import matplotlib.pyplot as plt
import numpy as np
import sys

h_minexp = int(sys.argv[1])
h_maxexp = int(sys.argv[2])
h_norms_file = sys.argv[3]
h_plot_file = sys.argv[4]

# Define the results directory
results_dir = "results"

h_vals = [10**h for h in range(h_minexp, h_maxexp + 1)]

data_matrix = np.loadtxt(h_norms_file)   # Shape: hvals x 3

# Plotting
fig, ax = plt.subplots(figsize=(5,3))

beta_labels = [r'$\beta^1$', r'$\beta^2$', r'$\beta^3$'] 

# Plot each row in data_matrix (norms for each beta term)
Norders = data_matrix.shape[1]
colors = ['tab:blue', 'tab:orange', 'tab:green']
for (i,color) in zip(range(0,Norders), colors):
    ax.plot(h_vals, data_matrix[:,i], marker='o', label=beta_labels[i], color=color)

# Set axes and labels
ax.set_xscale('log')  # Logarithmic scale for h values
ax.set_yscale('log')  # Logarithmic scale for norm values
ax.set_xlabel(r"h / Ha$^{-1}$")
ax.set_ylabel("Relative error")
ax.legend(title="Beta Terms", loc='best')

plt.tight_layout()

# Save the plot as a PNG file in the results directory
plt.savefig(h_plot_file, dpi=300, bbox_inches='tight')
