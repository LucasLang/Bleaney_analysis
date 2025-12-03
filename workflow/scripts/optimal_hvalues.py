import numpy as np
import sys

h_minexp = int(sys.argv[1])
h_maxexp = int(sys.argv[2])
h_norms_file = sys.argv[3]
output_file = sys.argv[4]

# Define the results directory
results_dir = "results"

h_vals = [10**h for h in range(h_minexp, h_maxexp + 1)]

data_matrix = np.loadtxt(h_norms_file)   # Shape: hvals x 3

# Assuming h_vals is initially a list
h2_vals = np.array(h_vals)  # Convert to NumPy array for indexing

# Find index of minimum for each beta term
min_indices = np.argmin(data_matrix, axis=0).flatten()  # Ensure indices are 1D

# Get h values at minimum error
optimal_h_values = h2_vals[min_indices.astype(int)]

# Compute the average of the optimal h values
optimal_h_avg = np.mean(optimal_h_values)

np.savetxt(output_file, [optimal_h_avg])

