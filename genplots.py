# from matplotlib import pyplot
import numpy as np
from brainnets import aux


def plot_mean_with_shaded_errors(ax, ydata, conf_int, xdata=None,
                                 xticklabels=None,
                                 color=None,
                                 label=None):
    """
    Parameters
    ----------
    ax : matplotlib axes object
        axes to be plotted to
    ydata : numpy array
        the mean values
    conf_int : tuple of numpy arrays (or lists)
        ``conf_int[0]`` contains the lower values
        ``conf_int[1]`` contains the higher values
    xdata : numpy array of values
    xticklabels : numpy array, list
        list of tick-values for the x-axis
    color : a matplotlib color
        what color to use for plotting
    label : str
        the label of the mean + error
    """

    if xdata is None:
        xdata = range(len(ydata))
    if label is not None:
        ax.plot(xdata, ydata, color=color, lw=1.5, label=label)
    else:
        ax.plot(xdata, ydata, color=color, lw=1.5)
    ax.set_xlim((np.min(xdata), np.max(xdata)))
    alpha = 0.25
    if xticklabels is None:
        xticklabels = ax.get_xticklabels()
    ax.fill_between(xdata, conf_int[0], conf_int[1], color=color, alpha=alpha)


def plot_inv_cdf(ax, data, xscale="linear", yscale="linear", **plottingArgs):
    """
    Plots the 1-CDF for data.

    Currently only a naive implementation available.
    (Bad implementation with a lot of data points)

    Parameters
    ----------
    ax : matplotlib axes object
        the axes to plot to
    data : 1D numpy array

    xscale : {"linear", "log"}, optional
        the xscale in the plot
    yscale : {"linear", "log"}, optional
        the xscale in the plot
    **plottingArgs : dict
        any keyword arguments that can be passed on to matplotlib.pyplot.plot
        function
    """
    # sort values
    data.sort()  # from small to big
    probs = np.linspace(1., 1. / len(data), len(data))
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.plot(data, probs, **plottingArgs)
    ax.legend()
    return ax


def plot_lin_pdf(ax, values, n_bins=100, color="r",
                 label=None, normalize=True):
    """
    Plots a PDF with linear bins.

    Parameters
    ----------
    ax : matplotlib axes object
        the axes to plot to
    values : 1D numpy array
        the values of the pdf
    n_bins : int, optional
        the number of bins in the pdf
    color : any valid matplotlib color, optional
        the color which the pdf is plotted with
    label : str, optional
        the label of the pdf
    normalize : bool, optional
        if False, the counts are shown instead of the pdf
    """
    values = np.array(values)
    val_max = np.max(values)
    val_min = np.min(values)
    linbins, bin_centers = aux.get_lin_bins(
        n_bins,
        val_min - 0.1 * np.abs(val_min),
        val_max + 0.1 * np.abs(val_max)
    )
    bin_counts = np.bincount(
        np.digitize(values, linbins), minlength=n_bins + 1)[1:]
    # normalize
    if normalize:
        bin_counts = bin_counts * 1. / \
            (len(values) * (linbins[1] - linbins[0]))
    if label is not None:
        ax.plot(bin_centers, bin_counts, color=color, label=label)
    else:
        ax.plot(bin_centers, bin_counts, color=color)
