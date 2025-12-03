import numpy as np
import os
import sys

h_minexp = int(sys.argv[1])
h_maxexp = int(sys.argv[2])

# Define the results directory
results_dir = "results"

h_vals = [10**h for h in range(h_minexp, h_maxexp + 1)]

data_matrix = np.loadtxt(os.path.join(results_dir, "h_norms.txt"))   # Shape: hvals x 3

# Assuming h_vals is initially a list
h2_vals = np.array(h_vals)  # Convert to NumPy array for indexing

# Find index of minimum for each beta term
min_indices = np.argmin(data_matrix, axis=0).flatten()  # Ensure indices are 1D

# Get h values at minimum error
optimal_h_values = h2_vals[min_indices.astype(int)]

# Compute the average of the optimal h values
optimal_h_avg = np.mean(optimal_h_values)

output_file = os.path.join(results_dir, "optimal_h_value.txt")
np.savetxt(output_file, [optimal_h_avg])
print(f"Results saved to {output_file}")

