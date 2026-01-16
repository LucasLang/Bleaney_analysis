import numpy as np
import sys
from py_functions import barplot_chi

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

methods = ["Bleaney (2nd order)", "3rd order", "Exact dyadic", "CASSCF"]
ylabel = r"$\chi_\mathrm{ax}$ / $\AA^3$"
barplot_chi(Ln_list, methods, values, ylabel, outfile)

