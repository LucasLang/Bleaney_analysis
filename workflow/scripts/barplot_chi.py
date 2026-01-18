import matplotlib.pyplot as plt
import numpy as np
import sys
from py_functions import barplot_chi

Ln_list_file = sys.argv[1]
chi_files = [sys.argv[2], sys.argv[3]]
chi_CASSCF_files = [sys.argv[4], sys.argv[5]]
outfile = sys.argv[6]

with open(Ln_list_file, "r") as f:
    Ln_list = [line.strip() for line in f]

values_calculated = [np.loadtxt(file) for file in chi_files]
values_CASSCF = [np.loadtxt(file) for file in chi_CASSCF_files]
values_CASSCF = [np.atleast_2d(v).T for v in values_CASSCF]
values = [np.hstack((values_calculated[0], values_CASSCF[0])),
          np.hstack((values_calculated[1], values_CASSCF[1]))]

methods = ["Bleaney (2nd order)", "3rd order", "Exact dyadic", "CASSCF"]
ylabel_ax = r"$\chi_\mathrm{ax}$ / $\AA^3$"
ylabel_rh = r"$\chi_\mathrm{rh}/\chi_\mathrm{ax}$"


fig = plt.figure(figsize=(6.7, 3), layout="constrained")
gs = fig.add_gridspec(nrows=2, ncols=2, height_ratios=[0.18, 1.0])

axes = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
ax_legend = fig.add_subplot(gs[0, :])
ax_legend.axis("off")


barplot_chi(axes[0], Ln_list, methods, values[0], ylabel_ax, legend=False)   # chi_ax
barplot_chi(axes[1], Ln_list, methods, values[1], ylabel_rh, legend=False)   # chi_rh

# Global legend
handles, labels = axes[0].get_legend_handles_labels()     # same labels in both plots
ax_legend.legend(handles, labels,
          loc='center',
          #bbox_to_anchor=(0.01, 0.09),
          borderaxespad=0,
          frameon=True,
          ncol = 4)

plt.tight_layout()
plt.savefig(outfile, dpi=300)
