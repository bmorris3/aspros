import numpy as np
from scipy.signal import find_peaks

__all__ = ['bls_peakfinder']


def bls_peakfinder(results):
    """
    Find peaks in a Box Least Squares spectrum.

    Parameters
    ----------
    results : `~astropy.timeseries.BoxLeastSquaresResults`
        BLS results

    Returns
    -------
    inds : `~numpy.ndarray`
        Indices with the top powers, sorted in order of peak height
    significance : float
        Ratio of the height of the tallest peak to the height of the
        second tallest peak
    """
    maxima = find_peaks(results.power, distance=100)[0]

    top_power_inds = maxima[np.argsort(results.power[maxima])[::-1]]

    highest_peak = results.power[top_power_inds[0]]
    next_highest_peak = results.power[top_power_inds[1]]

    significance = highest_peak / next_highest_peak

    return top_power_inds, significance