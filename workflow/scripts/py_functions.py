import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

def barplot_chi(rowlabels, columnlabels, data, ylabel, outfile):
    x = np.arange(len(rowlabels))  # this is a dummy number that is not displayed in the plot
    width = 0.18  # width of each bar

    fig, ax = plt.subplots(figsize=(5, 3))

    # Colors for each method
    colors = ['tab:orange', 'tab:green', 'tab:red', 'tab:purple']

    # Plot each method's bars
    for i, (method, color) in enumerate(zip(columnlabels, colors)):
        ax.bar(x + i * width - 1.5 * width, data[:, i], width, label=method, color=color, edgecolor='black', linewidth=0.8)

    # Axis and labels
    ax.axhline(0, color='black', linewidth=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(rowlabels)
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False)
    ax.set_xlim(x[0] - 0.5, x[-1] + 0.5)
    ax.tick_params(direction='in', top=True, right=True)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)

def plot_error_vs_h(ax, h_vals, h_norms_file, legend=True):
    data_matrix = np.loadtxt(h_norms_file)   # Shape: hvals x 3

    beta_labels = ['1', '2', '3'] 

    # Plot each row in data_matrix (norms for each beta term)
    Norders = data_matrix.shape[1]
    colors = ['tab:blue', 'tab:orange', 'tab:green']
    for (i,color) in zip(range(0,Norders), colors):
        ax.plot(h_vals, data_matrix[:,i], marker='o', label=beta_labels[i], color=color)

    # Set axes and labels
    ax.set_xscale('log')  # Logarithmic scale for h values
    ax.set_yscale('log')  # Logarithmic scale for norm values
    ax.set_xlabel(r"$h$ / Ha$^{-1}$")
    ax.set_ylabel("Relative error")
    ax.set_xticks(h_vals)
    if legend:
        ax.legend(title="Derivative order", loc='best')

def convergence_plot(ax, beta_terms, data_matrix, formatted_labels, temps, norm):
    cmap = cm.viridis
    colors = [cmap(norm(t)) for t in temps]

    # Plotting with improved visibility
    for i, (temp_data, color) in enumerate(zip(data_matrix, colors)):
        ax.plot(beta_terms, temp_data, marker='o', label=formatted_labels[i], 
                color=color)
    ax.set_yscale('log')
    ax.set_xlabel("Highest-order term in the expansion")
    ax.set_ylabel("Relative error")
    ax.grid(True, which="both", ls="-", alpha=0.2)

def sci_label(x, sig=2):
    """Return a MathText string like r'$m\times10^{e}$' for a number x."""
    if x == 0 or not np.isfinite(x):
        return r"$0$"
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / (10 ** exp)
    # Keep the sign with the mantissa
    return rf"${mant:.{sig}f}\times 10^{{{exp}}}$"
