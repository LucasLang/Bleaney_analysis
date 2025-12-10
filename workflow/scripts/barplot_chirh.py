import numpy as np
import matplotlib.pyplot as plt
import sys

Ln_list_file = sys.argv[1]
chi_ax_file = sys.argv[2]
chi_ax_CASSCF_file = sys.argv[3]
outfile = sys.argv[4]

with open(Ln_list_file, "r") as f:
    Ln_list = [line.strip() for line in f]

values_calculated = np.loadtxt(chi_ax_file)
values_CASSCF = np.loadtxt(chi_ax_CASSCF_file)
values_CASSCF = np.atleast_2d(values_CASSCF).T
values = np.hstack((values_calculated, values_CASSCF))

methods = ["Bleaney", "3rd order", "exact", "CASSCF"]

x = np.arange(len(Ln_list))  # this is a dummy number that is not displayed in the plot
width = 0.18  # width of each bar

fig, ax = plt.subplots(figsize=(8, 5))

# Colors for each method
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

# Plot each method's bars
for i, (method, color) in enumerate(zip(methods, colors)):
    ax.bar(x + i * width - 1.5 * width, values[:, i], width, label=method, color=color, edgecolor='black', linewidth=0.8)

# Axis and labels
ax.axhline(0, color='black', linewidth=0.9)
ax.set_xticks(x)
ax.set_xticklabels(Ln_list)
ax.set_ylabel("Calculated Property (units)")
ax.set_xlabel("Molecule")
ax.legend(frameon=False)
ax.set_xlim(x[0] - 0.5, x[-1] + 0.5)
ax.tick_params(direction='in', top=True, right=True)

plt.tight_layout()
plt.savefig(outfile, dpi=300)
