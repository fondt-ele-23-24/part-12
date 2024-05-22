from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
from matplotlib import pyplot as plt
from matplotlib import cm

def plot(x, y, xlabel=None, ylabel=None, title=None, figsize=None):
    plt.figure(figsize=figsize)
    plt.plot(x, y)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=14)
    plt.grid(':')
    plt.show()

    
def plot_univariate_function(f, x, figsize=None):
    plt.figure(figsize=figsize)
    plt.plot(x, f(x))
    plt.plot(plt.xlim(), [0, 0])
    plt.grid()
    plt.show()


def plot_state_evolution(X, t, ylabels=None, xlabel=None, title=None, figsize=None, same_scale=False):
    # Costruisco una figura con diverse sotto-figure
    fig = plt.figure(figsize=figsize)

    # Preparo la mappa dei colori
    cmap = cm.get_cmap('Set2')

    # Define font size
    fontsize = None if figsize is None else 0.9 * figsize[0]

    # Considero le colonne con le diverse componenti dello stato
    n = X.shape[1]
    for i in range(0, n):
        if not same_scale:
            plt.subplot(n, 1, i+1)
        # Draw a curve
        y_lbl = ylabels[i] if ylabels is not None else None
        plt.plot(t, X[:, i], color=cmap(i), label=y_lbl)
        # Add y label if not using the same scale
        if not same_scale:
            plt.ylabel(y_lbl)
        # Label the x axis
        plt.xlabel(xlabel)
        # Add a grid
        plt.grid(':')

    # Define the title
    if same_scale:
        plt.title(title)
    else:
        plt.suptitle(title)

    # Crop useless space on the borders
    plt.tight_layout()

    if same_scale:
        plt.legend()
        
    plt.show()
