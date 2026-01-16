import matplotlib.pyplot as plt
import numpy as np

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
