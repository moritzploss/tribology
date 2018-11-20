# -*- coding: utf-8 -*-

"""

This module contains functions related to Hertz contact theory. As of now, the
equations are limited to elliptical and circular contacts. Equations for line
contacts (cylinder-on-flat or cylinder-on-cylinder) are currently not
implemented.

"""
import copy
import warnings
from math import sqrt, pi, log, floor

import numpy as np

from tribology.boundary_element import __secant
from tribology.tribology import profball


def __rred(r_1, r_2):
    """

    Calculate the reduced (effective) radius of two radii according to Hertzian
    contact theory.

    Parameters
    ----------
    r_1: scalar
        The first radius.
    r_2: scalar
        The second radius.

    Returns
    -------
    r_red: scalar
        The reduced (effective) radius.

    """
    if (r_1 == 0 or abs(r_1) == float('inf')) and r_2 != 0:
        r_red = r_2
    elif (r_2 == 0 or abs(r_2) == float('inf')) and r_1 != 0:
        r_red = r_1
    elif (r_1 == 0 or abs(r_1) == float('inf')) and \
            (r_2 == 0 or abs(r_2) == float('inf')):
        r_red = 0
    elif r_1 == -r_2:
        r_red = 0
    else:
        r_red = 1 / (1 / r_1 + 1 / r_2)
    return r_red


def reff(r_x_1, r_y_1, r_x_2, r_y_2, ):
    """

    Calculate the effective radii for two bodies according to Hertzian contact
    theory. It is assumed that the two major axis of each body (x- and y-axis)
    are perpendicular to each other and that the x and y axes of both bodies are
    aligned.

    Parameters
    ----------
    r_x_1: scalar
        The radius of body 1 in direction 1 (x).
    r_y_1: scalar
        The radius of body 1 in direction 2 (y).
    r_x_2: scalar
        The radius of body 2 in direction 1 (x).
    r_y_2: scalar
        The radius of body 2 in direction 2 (y).

    Returns
    -------
    r_eff: scalar
        The effective radius.
    r_eff_x: scalar
        The effective radius in x-direction.
    r_eff_y: scalar
        The effective radius in y-direction.

    """
    recip_radius = []
    for radius in [r_x_1, r_y_1, r_x_2, r_y_2]:
        if radius == 0:
            recip_radius.append(float('inf'))
        else:
            recip_radius.append(1 / radius)

    r_eff_x = __rred(r_x_1, r_x_2)
    r_eff_y = __rred(r_y_1, r_y_2)
    r_eff = __rred(r_eff_x, r_eff_y)
    return r_eff, r_eff_x, r_eff_y


def eeff(e_1, nu_1, e_2, nu_2):
    """

    Calculate the effective (Young's) modulus of two contact bodies according
    to Hertzian contact theory.

    Parameters
    ----------
    e_1: ndarray, scalar
        The Young's modulus of contact body 1.
    nu_1: ndarray, scalar
        The Poisson ratio of contact body 1.
    e_2: ndarray, scalar
        The Young's modulus of contact body 2.
    nu_2: ndarray, scalar
        The Poisson ratio of contact body 1.

    Returns
    -------
    e_eff: scalar
        The effective modulus.

    """
    e_eff = 1 / ((1 - nu_1 ** 2) / (2 * e_1) + (1 - nu_2 ** 2) / (2 * e_2))
    return e_eff


def __auxparamshertz(r_eff_x, r_eff_y):
    """

    Calculate a set of parameters required for Hertz contact area/pressure
    calculations.

    Parameters
    ----------
    r_eff_x: scalar
        The effective redius of the contact problem in x-direction.
    r_eff_y: scalar
        The effective redius of the contact problem in x-direction.

    Returns
    -------
    a_ast: scalar
        A parameter often referred to as a*.
    b_ast: scalar
        A parameter often referred to as b*.
    kappa: scalar
        A parameter often referred to as kappa.
    param_lambda: scalar
        A parameter often referred to as lambda.

    """
    if r_eff_x == 0 or r_eff_y == 0:
        param_lambda = 0
    else:
        param_lambda = min(r_eff_x / r_eff_y, r_eff_y / r_eff_x)

    if param_lambda == 0:
        param_lambda = 10**-256

    kappa = 1 / (1 + sqrt(log(16 / param_lambda) / (2 * param_lambda)) -
                 sqrt(log(4)) + 0.16 * log(param_lambda))
    a_ast = kappa * (1 + (2 * (1 - kappa ** 2)) / (pi * kappa ** 2) -
                     0.25 * log(kappa)) ** (1 / 3)
    b_ast = a_ast / kappa
    return a_ast, b_ast, kappa, param_lambda


def dhertz(e_eff, r_x_1, r_y_1, r_x_2, r_y_2, force):
    """

    Calculate the elastic normal displacement in an elliptic (including
    circular) Hertzian contact.

    Parameters
    ----------
    e_eff: scalar
        The effective modulus of the contact problem.
    r_x_1: scalar
        The radius of body 1 in direction 1 (x).
    r_y_1: scalar
        The radius of body 1 in direction 2 (y).
    r_x_2: scalar
        The radius of body 2 in direction 1 (x).
    r_y_2: scalar
        The radius of body 2 in direction 2 (y).
    force: scalar
        The normal force in the contact.

    Returns
    -------
    norm_disp: scalar
        The elastic normal displacement of the contact problem.

    """
    apb = 0.5 * (1 / r_x_1 + 1 / r_y_1 + 1 / r_x_2 + 1 / r_y_2)
    bma = 0.5 * sqrt((1 / r_x_1 - 1 / r_y_1) ** 2 +
                     (1 / r_x_2 - 1 / r_y_2) ** 2)
    r_a = 1 / (apb - bma)
    r_b = 1 / (apb + bma)
    r_c = sqrt(r_a * r_b)

    f_1 = 1 - pow(pow(r_a / r_b, 0.0602) - 1, 1.456)
    f_2 = 1 - pow(pow(r_a / r_b, 0.0684) - 1, 1.531)

    param_c = pow(3 * force * r_c / (4 * e_eff), 1 / 3) * f_1
    param_e = 1 - pow(r_b / r_a, 4 / 3)
    param_a = param_c * pow(1 - param_e ** 2, 1 / 4)
    param_b = param_c * pow(1 - param_e ** 2, 1 / 4)
    norm_disp = param_a * param_b / r_c * (f_2 / f_1)
    return norm_disp


def ahertz(r_eff, r_eff_x, r_eff_y, e_eff, force):
    """

    Calculate the contact area in an elliptic (including circular) Hertzian
    contact.

    Parameters
    ----------
    r_eff: scalar
        The effective radius of the contact problem.
    r_eff_x: scalar
        The effective radius of the contact problem in x-direction.
    r_eff_y: scalar
        The effective radius of the contact problem in y-direction.
    e_eff: scalar
        The effective modulus of the contact problem.
    force: scalar
        The normal force in the contact.

    Returns
    -------
    half_axis_a: scalar
        The contact area half axis in x-direction.
    half_axis_b: scalar
        The contact area half axis in y-direction.
    a_hertz: scalar
        The contact area.

    """
    a_ast, b_ast, _, _ = __auxparamshertz(r_eff_x, r_eff_y)
    half_axis_a = a_ast * (3 * force * r_eff / e_eff) ** (1 / 3)
    half_axis_b = b_ast * (3 * force * r_eff / e_eff) ** (1 / 3)
    a_hertz = pi * half_axis_a * half_axis_b
    return half_axis_a, half_axis_b, a_hertz


def fhertz(r_eff, r_eff_x, r_eff_y, e_eff, p_critical):
    """

    Calculate the load carrying capacity of an elliptic (including circular)
    Hertzian contact.

    Parameters
    ----------
    r_eff: scalar
        The effective radius of the contact problem.
    r_eff_x: scalar
        The effective radius of the contact problem in x-direction.
    r_eff_y: scalar
        The effective radius of the contact problem in y-direction.
    e_eff: scalar
        The effective modulus of the contact problem.
    p_critical: scalar
        The critical mean Hertzian contact pressure that the contact can sustain
        without plastic deformation.

    Returns
    -------
    f_crit: scalar
        The load carrying capacity of the Hertzian contact.

    """
    a_ast, b_ast, _, _ = __auxparamshertz(r_eff_x, r_eff_y)
    f_crit = (pi * a_ast * b_ast * p_critical) ** 3 * (3 * r_eff / e_eff) ** 2
    return f_crit


def phertz(r_eff, r_eff_x, r_eff_y, e_eff, force, ret='mean'):
    """

    Calculate the contact pressure in an elliptic (including circular) Hertzian
    contact. By default, the mean pressure is returned. Use :code:`ret='max'`
    for maximum pressure.

    Parameters
    ----------
    r_eff: scalar
        The effective radius of the contact problem.
    r_eff_x: scalar
        The effective radius of the contact problem in x-direction.
    r_eff_y: scalar
        The effective radius of the contact problem in y-direction.
    e_eff: scalar
        The effective modulus of the contact problem.
    force: scalar
        The normal force in the contact.
    ret: str, optional
        Return mean Hertzian pressure for :code:`'mean'` and maximum pressure
        for :code:`'max'`

    Returns
    -------
    p_hertz: scalar
        The Hertzian contact pressure in the contact.

    """
    _, _, area = ahertz(r_eff, r_eff_x, r_eff_y, e_eff, force)
    p_hertz = force / area

    if ret == 'mean':
        pass
    elif ret == 'max':
        p_hertz *= 1.5
    else:
        raise ValueError("argument 'ret' must have value 'max' or 'mean'")

    return p_hertz


def __circ3points(x_cords, y_cords):
    """

    Calculate the radius of a circle that goes through three arbitratry
    points in the x-y plane.

    Parameters
    ----------
    x_cords: tuple
        x-coordinates of the three points: (x1, x2, x3)
    y_cords: tuple
        y-coordinates of the three points: (y1, y2, y3)

    Returns
    -------
    rad: float
        The radius of the circle. If all points are on a straight line, rad
        will be :code:`inf`

    """
    a = sqrt((x_cords[0] - x_cords[1]) ** 2 + (y_cords[0] - y_cords[1]) ** 2)
    b = sqrt((x_cords[1] - x_cords[2]) ** 2 + (y_cords[1] - y_cords[2]) ** 2)
    c = sqrt((x_cords[2] - x_cords[0]) ** 2 + (y_cords[2] - y_cords[0]) ** 2)
    s = (a + b + c) / 2
    d = sqrt(s * (s - a) * (s - b) * (s - c))

    if d != 0:
        rad = a * b * c / (4 * d)
    else:
        rad = float('Inf')

    return rad


def approx_hertz_rad(axis, profile, iterations=10):
    """

    Approximate the Hertz contact radius of a profile, i.e., find the radius
    of a circle that minimizes the average discrepancy between the actual
    profile and a circle profile. The secant method is used to find the circle
    approximation.

    This method has been tested for continuous (or at least semi-continuous)
    profiles only. It works best for cases where the circle radius is larger
    than the absolute minimum and maximum of the profile axis.

    The function will not produce sensible results for profiles where the
    maximum profile height is large compared to the axis length.

    Parameters
    ----------
    axis: ndarray
        The coordinates of the `profile` points along their reference axis. For
        example, if the `profile` array contains the profile heights along the
        x-axis in the range [-5, 5], then `axis` contains the corresponding
        x-axis values.
    profile: ndarray
        The profile heights along `axis`.
    iterations: int, optional
        Number of times the secant method is applied before the radius is
        returned.

    Returns
    -------
    rad: scalar
        The radius of the circle that best approximates `profile`. If the
        profile cannot be approximated with a circle (usually if the profile
        is a straight line), the function returns :code:`inf`.
    ax_approx: ndarray
        The elements of :code:`axis` that were used for the approximation.
    prof_approx: ndarray
        The elements of :code:`profile` that were used for the approximation.

    """

    if abs(np.amax(profile) - np.amin(profile)) > \
            abs(np.amax(axis) - np.amin(axis)):
        warnings.warn('profile range larger than axis range')

    cen_idx = floor(len(axis) / 2)

    # flip profile so that lowest point is not at the edges and minimum at zero
    if (profile[cen_idx] > profile[0] and profile[cen_idx] > profile[-1]) or \
            profile[0] > profile[cen_idx] > profile[-1]:
        profile *= -1
    profile -= np.amin(profile)

    # approximate profile with circle for initial guess of rad
    if profile[0] <= profile[-1]:
        idx = np.argmin(np.abs(profile[cen_idx:-1] - profile[0])) - 1
    else:
        idx = np.argmin(np.abs(profile[0:cen_idx] - profile[-1])) + 1
    rad = __circ3points(
        (axis[0], axis[cen_idx], axis[idx]),
        (profile[0], profile[cen_idx], profile[idx]))

    ax_approx = copy.deepcopy(axis)
    prof_approx = copy.deepcopy(profile)

    # improve guess with secant method
    if rad != float('Inf'):
        radii = []
        deltas = []
        for _ in range(iterations):
            ax_approx = copy.deepcopy(axis)
            prof_approx = copy.deepcopy(profile)

            # remove axis elements that lie outside the circle
            while abs(ax_approx[-1]) >= rad:
                idx = len(ax_approx) - 1
                ax_approx = np.delete(ax_approx, idx)
                prof_approx = np.delete(prof_approx, idx)
            while abs(prof_approx[0]) >= rad:
                ax_approx = np.delete(ax_approx, 0)
                prof_approx = np.delete(prof_approx, 0)

            # apply secant method
            radii = np.append(radii, [rad])
            delta = np.sum(profball(ax_approx, rad) - prof_approx)
            deltas = np.append(deltas, [delta])
            rad = __secant(radii, deltas)
            if rad < 0:
                rad *= -1

    return rad, ax_approx, prof_approx
