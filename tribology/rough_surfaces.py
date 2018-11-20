# -*- coding: utf-8 -*-

"""

This module contains functions related to rough surfaces. Please see the
:code:`/demo` directory for minimum examples.

"""


import numpy as np


def abbottfirestone(trace, num_bins=100):
    """

    Calculate the Abbott-Firestone curve for a 1D profile trace.

    Parameters
    ----------
    trace: ndarray
        The profile heights of the 1D trace.
    num_bins: positive int, optional
        The number of bins for the calculation of the profile height
        probability distribution.

    Returns
    -------
    bins: ndarray
        The bins used for the calculation of the profile height
        probability distribution (= x-axis data for Abbott-Firestone plot).
    prob_dist: ndarray
        The probability distribution of the profile (= y-axis data for Abbott-
        Firestone plot).

    """
    bins = np.linspace(np.amax(trace), np.amin(trace), num_bins)
    prob_dist = []
    for each_bin in bins:
        prob_dist.append(
            len(np.where(trace >= each_bin)[0]) / len(trace) * num_bins)
    return bins, prob_dist


def randsurf(n_x, n_y, delta_x, delta_y, s_q, lambda_x, lambda_y):
    """

    Generate a random rough surface in the x-y-plane. The surface has a Gaussian
    roughness distribution, a mean surface roughness of (approximately) 0, and
    a standard deviation of :code:`s_a`.

    The correlation length can be specified separately in x and y; however,
    the correlation length and grid spacing need to be constant in both
    directions. The computational cost increases with increasing correlation
    length in either direction.

    The implementation follows the approach described in the appendix of:

    N. Patir and H. S. Cheng, An Average Flow Model for Determining Effects
    of Three-Dimensional Roughness on Partial Hydrodynamic Lubrication, Journal
    of Lubrication Technology, 100(1), pages 12-17, 1978. doi:10.1115/1.3453103

    The implementation is similar to that used by Björklund and Andersson:

    S. Björklund and S. Andersson, A numerical method for real elastic contacts
    subjected to normal and tangential loading, Wear, 179, pages 117-122, 1994.
    doi: 10.1016/0043-1648(94)90228-3

    Parameters
    ----------
    n_x: int
        Number of grid points along x-axis.
    n_y: int
        Number of grid points along y-axis.
    delta_x: float
        Distance between grid points along the x-axis
    delta_y: float
        Distance between grid points along the y-axis.
    s_q: float
        Standard deviation of the surface roughness.
    lambda_x: float
        Correlation length along the x-axis.
    lambda_y: float
        Correlation length along the y-axis.

    Returns
    -------
    heights: ndarray
        Surface heights in the x-y-plane.

    """
    m = int(np.round(2 * lambda_x / delta_x))
    n = int(np.round(2 * lambda_y / delta_y))

    if m < 1 or n < 1:
        raise ValueError('ratio between correlation length and grid spacing '
                         'too small. increase lambda or decrease delta.')
    else:
        surf_rand = np.random.randn(n_x + m, n_y + n)
        heights = np.zeros((n_x, n_y))

        for i in range(0, n_x):
            for j in range(0, n_y):
                heights[i, j] = np.sum(
                    surf_rand[i + 1:i + m + 1, j + 1:j + n + 1])

        # set std according to s_q argument, and mean equal to zero
        heights = heights / np.std(heights) * s_q
        heights -= np.mean(heights)

    return heights
