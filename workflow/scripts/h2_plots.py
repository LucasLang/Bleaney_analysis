import matplotlib.pyplot as plt
import numpy as np
import json
import os

# Define the results directory
results_dir = "results"

# Read h values from hruns.json (in the results directory)
hruns_file = os.path.join(results_dir, "hruns.json")
with open(hruns_file, "r") as file:
    hruns = json.load(file)

h_vals = [h[1] for h in hruns]  # Extract h values 

# Define file paths in the 'results' folder
file_paths = [
    os.path.join(results_dir, "h_norm1.txt"),
    os.path.join(results_dir, "h_norm2.txt"),
    os.path.join(results_dir, "h_norm3.txt")
    # os.path.join(results_dir, "h_norm4.txt"),
    # os.path.join(results_dir, "h_norm5.txt"),
    # os.path.join(results_dir, "h_norm6.txt")
]

# Load data into a matrix
data_matrix = []
for path in file_paths:
    data = np.loadtxt(path)
    data_matrix.append(data)

data_matrix = np.array(data_matrix)  # Shape: (6, len(h_vals))

# Assuming h_vals is initially a list
h2_vals = np.array(h_vals)  # Convert to NumPy array for indexing

# Find index of minimum for each beta term
min_indices = np.argmin(data_matrix, axis=1).flatten()  # Ensure indices are 1D

# Get h values at minimum error
optimal_h_values = h2_vals[min_indices.astype(int)]

# Compute the average of the optimal h values
optimal_h_avg = np.mean(optimal_h_values)

# Save results to a text file
output_file = os.path.join(results_dir, "optimal_h_values.txt")
with open(output_file, "w") as f:
    f.write("Optimal h values for individual beta terms:\n")
    for i, h in enumerate(optimal_h_values):
        f.write(f"Beta term {i+1}: {h:.15e}\n")
    f.write(f"\nOptimal h value (average of all beta terms): {optimal_h_avg:.15e}\n")

print(f"Results saved to {output_file}")


# Plotting
fig, ax = plt.subplots(figsize=(10, 6))

beta_labels = [r'$\beta^1$', r'$\beta^2$', r'$\beta^3$'] #, r'$\beta^4$', r'$\beta^5$', r'$\beta^6$'

# Plot each row in data_matrix (norms for each beta term)
for i, norm_values in enumerate(data_matrix):  # Rows correspond to beta terms
    ax.plot(h_vals, norm_values, marker='o', label=beta_labels[i])

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
np.savetxt(plot_data_file, np.column_stack((h_vals, data_matrix.T)), header='h_vals beta1 beta2 beta3')

# Save the plot as a PNG file in the results directory
h_plot_file = os.path.join(results_dir, "h_plot.png")
plt.savefig(h_plot_file, dpi=300, bbox_inches='tight')