import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import sys
from fractions import Fraction

def sci_label(x, sig=2):
    """Return a MathText string like r'$m\times10^{e}$' for a number x."""
    if x == 0 or not np.isfinite(x):
        return r"$0$"
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / (10 ** exp)
    # Keep the sign with the mantissa
    return rf"${mant:.{sig}f}\times 10^{{{exp}}}$"

k = 3.166811563e-6 
conversion_factor = 219474.63
kB = k * conversion_factor

S = float(sys.argv[1])
E = int(sys.argv[2])
D = int(sys.argv[3])
diff_norms_file = sys.argv[4]
Tpointsfile = sys.argv[5]
error_analysis_plot_file = sys.argv[6]


temps = np.loadtxt(Tpointsfile)

data_matrix = np.loadtxt(diff_norms_file)

# X-axis labels
beta_terms = [r'$\beta$', r'$\beta^2$', r'$\beta^3$', r'$\beta^4$', r'$\beta^5$', r'$\beta^6$']

# Create figure and axes with larger size
plt.figure(figsize=(10,6))

ax = plt.gca()

# --- colors by temperature (log-normalized) ---
norm = mcolors.LogNorm(vmin=min(temps), vmax=max(temps))
cmap = cm.viridis
colors = [cmap(norm(t)) for t in temps]

# Plotting with improved visibility
for i, (temp_data, color) in enumerate(zip(data_matrix, colors)):
    temperature_ratio = D / (kB * temps[i])
    formatted_label = sci_label(float(temperature_ratio), sig=2)
    ax.plot(beta_terms, temp_data, marker='o', label=formatted_label, 
            color=color, linewidth=2, markersize=8)

# Convert to fraction
E_over_D_fraction = Fraction(E,D)

twotimesS = int(2*S)
S_frac = Fraction(twotimesS, 2)  # Fraction objects nicely interpolate into strings

# Main plot formatting
ax.set_yscale('log')
ax.set_ylim(1e-14, 1e6)
ax.set_xlabel("Highest-order Term in the Expansion", fontsize=14, fontweight='bold')
ax.set_ylabel("Relative Error Values", fontsize=14, fontweight='bold')
#ax.set_title(f"TM", fontsize=14, fontweight='bold', pad=20)
ax.grid(True, which="both", ls="-", alpha=0.2)
ax.tick_params(axis='both', which='major', labelsize=15)

# First legend (D/kT values)
handles, labels = ax.get_legend_handles_labels()
first_legend = ax.legend(handles, labels, title="D/kT", 
          title_fontsize=14, fontsize=12,
          loc='lower left',
          bbox_to_anchor=(0.01, 0.09),
          borderaxespad=0,
          frameon=True,
          fancybox=True,
          shadow=True)



# Create dummy lines without actual lines for the second legend
dummy1 = plt.Line2D([], [], linestyle='none', label=f"E/D = {E_over_D_fraction}")
dummy2 = plt.Line2D([], [], linestyle='none', label=f"S = {S_frac}")

# Second legend (E/D as fraction and S)
secondary_legend = ax.legend(
    handles=[dummy1, dummy2],
    loc='lower left', 
    frameon=True, 
    fontsize=13
)

ax.add_artist(first_legend)  

plt.tight_layout()

# Save with higher resolution and larger size
plt.savefig(error_analysis_plot_file, dpi=300, bbox_inches='tight')