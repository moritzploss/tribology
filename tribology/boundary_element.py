# -*- coding: utf-8 -*-

"""

This module contains functions related to boundary element solvers for contact
mechanics calculations. Please see the :code:`/demo` directory for minimum
examples of how to combine the below functions to perform contact mechanics
calculations.

"""


import copy
from math import pi

import numexpr as ne
import numpy as np
import scipy.sparse.linalg as spla


def __secant(x_list, fx_list):
    """

    Applies secant method to find root *x* of function *f(x)*.
    If not enough (*x*, *f(x)*) value pairs are available to apply
    secant method, a new *x* value is guessed by slightly changing the
    initial *x* value

    Parameters
    ----------
    x_list: ndarray, list
        A list of known function inputs (x-values).
    fx_list: ndarray, list
        A list of known function outputs corresponding to the function inputs
        in ``x_list``.

    Returns
    -------
    x: scalar
        An estimate of the root of function f(x).

    """
    if fx_list[-1] != 0:
        if len(x_list) > 1 and \
                        fx_list[-1] != 0 and \
                        abs(fx_list[-1]) != abs(fx_list[-2]):
            x_0 = x_list[-2]
            x_1 = x_list[-1]
            fx_0 = fx_list[-2]
            fx_1 = fx_list[-1]
            slope = (fx_1 - fx_0) / (x_1 - x_0)
            return x_1 + (-fx_1 / slope)
        else:
            return x_list[0] * 0.9 + 0.0001
    else:
        return x_list[-1]


def __beinflumatgrid(axis):
    """

    Private helper method for influence matrix generation. Generates a grid
    that is required for the generation of an influence matrix.

    Parameters
    ----------
    axis: ndarray
        An axis of the elastic half space.

    Returns
    -------
    infl_mat_grid: ndarray
        A grid for influence matrix generation.

    """
    len_axis = len(axis)
    vec = np.zeros((1, len_axis))
    vec[0, :] = axis
    vertical_ax = np.zeros((len_axis, 1))
    vertical_ax[:, 0] = axis
    grid = np.repeat(vec, len_axis, axis=0)
    infl_mat_grid = np.absolute(np.subtract(grid, vertical_ax))
    return infl_mat_grid


def beinflumatred(infl_mat):
    """

    Calculate a reduced influence coefficient matrix from a complete influence
    coefficient matrix.

    Parameters
    ----------
    infl_mat: ndarray
        The complete influence coefficient matrix.

    Returns
    -------
    reduced_infl_mat: ndarray
        The reduced influence coefficient matrix. The matrix is square and of
        order ``n = np.shape(infl_mat)[0]`` :math:`\\times`
        ``np.shape(infl_mat)[1]``.

    """
    shape_mat = np.shape(infl_mat)
    len_mat = shape_mat[0] * shape_mat[1]
    reduced_infl_mat = np.zeros((len_mat, len_mat))
    counter = 0
    for i in range(0, shape_mat[0]):
        for j in range(0, shape_mat[1]):
            reduced_infl_mat[counter, :] = \
                np.reshape(infl_mat[i, j, :, :], len_mat)
            counter += 1
    return reduced_infl_mat


def beinflumat(x_axis, y_axis, e_eff):
    """

    Generate a complete influence coefficient matrix as required for boundary
    element contact mechanics calculations.

    Parameters
    ----------
    x_axis: ndarray
        The x-axis values of the coordinate grid.
    y_axis: ndarray
        The y-axis values of the coordinate grid.
    e_eff: scalar
        The effective modulus of the contact bodies.

    Returns
    -------
    infl_mat: ndarray
        The complete four-dimensional influence coefficient matrix of size
        ``len(x_axis)`` :math:`\\times` ``len(y_axis)`` :math:`\\times`
        ``len(x_axis)`` :math:`\\times` ``len(y_axis)``.

    """
    len_x = len(x_axis)
    len_y = len(y_axis)
    influence_matrix_complete = np.zeros((len_x, len_y, len_x, len_y))

    # generate coordinate grids
    a_factor = (x_axis[-1] - x_axis[0]) / (len_x - 1) / 2
    b_factor = (y_axis[-1] - y_axis[0]) / (len_y - 1) / 2
    x_grid = __beinflumatgrid(x_axis)
    y_grid = __beinflumatgrid(y_axis)

    # use numexpr to evaluate expressions
    xpa = ne.evaluate('x_grid + a_factor')
    xma = ne.evaluate('x_grid - a_factor')
    ypb = ne.evaluate('y_grid + b_factor')
    ymb = ne.evaluate('y_grid - b_factor')

    # calculate complete influence matrix
    for j in range(0, len_y):
        for j_prime in range(0, len_y):
            influence_matrix_complete[:, j, :, j_prime] =  \
                    (np.multiply(xpa, np.log(
                        np.divide(
                            ((ypb[j, j_prime]) +
                             np.sqrt(np.multiply((ypb[j, j_prime]),
                                                 (ypb[j, j_prime])) +
                                     np.multiply(xpa, xpa))),
                            ((ymb[j, j_prime]) +
                             np.sqrt(np.multiply((ymb[j, j_prime]),
                                                 (ymb[j, j_prime])) +
                                     np.multiply(xpa, xpa)))))) +
                     (ypb[j, j_prime]) * np.log(
                         np.divide(
                             (xpa +
                              np.sqrt(np.multiply((ypb[j, j_prime]),
                                                  (ypb[j, j_prime])) +
                                      np.multiply(xpa, xpa))),
                             (xma +
                              np.sqrt(np.multiply((ypb[j, j_prime]),
                                                  (ypb[j, j_prime])) +
                                      np.multiply(xma, xma))))) +
                     np.multiply(xma, np.log(
                         np.divide(
                             ((ymb[j, j_prime]) +
                              np.sqrt(np.multiply((ymb[j, j_prime]),
                                                  (ymb[j, j_prime])) +
                                      np.multiply(xma, xma))),
                             ((ypb[j, j_prime]) +
                              np.sqrt(np.multiply((ypb[j, j_prime]),
                                                  (ypb[j, j_prime])) +
                                      np.multiply(xma, xma)))))) +
                     (ymb[j, j_prime]) * np.log(
                         np.divide(
                             (xma +
                              np.sqrt(np.multiply((ymb[j, j_prime]),
                                                  (ymb[j, j_prime])) +
                                      np.multiply(xma, xma))),
                             (xpa +
                              np.sqrt(np.multiply((ymb[j, j_prime]),
                                                  (ymb[j, j_prime])) +
                                      np.multiply(xpa, xpa))))))

    infl_mat = influence_matrix_complete * 1 / e_eff * 2 / pi
    return infl_mat


def __begetd(profile, norm_disp):
    """

    Parameters
    ----------
    profile: ndarray
        The combined profiles of two contact bodies.
    norm_disp: scalar
        The normal elastic displacement (interference) of the two contact
        bodies.

    Returns
    -------
    disp: ndarray
        A displacement field containing the combined local elastic displacements
        of both contact bodies.

    """
    displ_field = np.subtract(np.ones(profile.shape) * norm_disp, profile)
    disp = np.reshape(displ_field, (profile.shape[0] * profile.shape[1]))
    return disp


def __solvepress(red_infl_mat, disp):
    """

    Solve for pressure distribution by removing negative pressure elements from
    the pressure matrix.

    Parameters
    ----------
    red_infl_mat: ndarray
        The recued influence matrix of the contact problem.
    disp: ndarray
        The (combined) displacement field of the contact problem.

    Returns
    -------
    pressure: ndarray
        The pressure field of the contact problem.

    """

    # find negative pressure arguments
    pressure = spla.gmres(red_infl_mat, disp)[0]
    p_index = np.zeros(len(pressure))
    negative_p = np.where(pressure < 0)[0]
    p_neg = copy.deepcopy(negative_p)

    while len(negative_p) > 0:
        pressure[p_neg] = 0
        p_index[p_neg] = 1
        u_new_reduced = np.delete(disp, [p_neg], axis=0)
        g_new_reduced = np.delete(red_infl_mat, [p_neg], axis=0)
        g_new_reduced = np.delete(g_new_reduced, [p_neg], axis=1)
        if pressure[np.where(p_index == 0)].size > 0:
            pressure[np.where(p_index == 0)] = \
                spla.gmres(g_new_reduced, u_new_reduced)[0]
        negative_p = np.where(pressure < 0)[0]
        p_neg = np.append(p_neg, negative_p)

    return pressure


def besolve(profile_1, profile_2, outer_force,
            red_infl_mat, delta_x, delta_y,
            norm_disp=0.1, max_offset=0.005):
    """

    Solve a system of linear equations to find the pressure and displacement
    distribution in a boundary element contact problem:

        :math:`[\\text{pressure}] = [\\text{reduced influence matrix}]^{-1}
        \\cdot [\\text{displacement}]`

    The start value for the relative normal displacement is defined by the
    ``norm_disp`` parameter. The default value corresponds to 0.1 units of
    length. If the start value of the normal displacement is close to the
    equilibrium displacement, the calculation time may be reduced significantly.

    The solver stops once inner and outer forces are in equilibrium. The maximum
    allowed discrepancy between inner and outer forces is defined by the
    ``max_offset`` parameter. The default value corresponds to 0.5 %.

    Parameters
    ----------
    profile_1: ndarray
        A matrix containing the profile heights of body 1.
    profile_2: ndarray
        A matrix containing the profile heights of body 2. The matrix must have
        the same size and grid spacing (in unit length) as ``profile_1``.
    outer_force: scalar
        The (outer) normal force acting on the contact.
    red_infl_mat: ndarray
        The reduced influence matrix calculated using the ``beinflumatred``
        method.
    delta_x: scalar
        The grid spacing of the profile matrices in unit length in x-direction.
    delta_y: scalar
        The grid spacing of the profile matrices in unit length in y-direction.
    norm_disp: scalar
        The initial normal elastic deformation in unit length used to start the
        solver.
    max_offset: scalar
        The maximum allowed percentage difference of inner and outer force (used
        to determine when solver has completed).

    Returns
    -------
    pressure: ndarray
        The contact pressure field in the contact.
    disp: ndarray
        The combined normal displacement field for contact body 1 and 2.
    inner_force: scalar
        The sum over all inner forces.
    norm_disp: scalar
        The normal elastic equilibrium displacement.

    """
    # initialise variables
    x_value = [0]
    fx_value = [outer_force]
    pressure = 0
    inner_force = 0
    profile = profile_1 + profile_2

    # while difference between inner forces and outer forces is significant
    while abs(inner_force - outer_force) > max_offset * outer_force:

        # update local displacements
        disp = __begetd(profile, norm_disp)
        pressure = __solvepress(red_infl_mat, disp)

        # calculate resulting force and adjust displacement for next loop
        pressure = np.reshape(pressure, (profile.shape[0], profile.shape[1]))
        inner_force = sum(sum(np.multiply(delta_x * delta_y, pressure)))
        x_value = np.append(x_value, [norm_disp])
        fx_value = np.append(fx_value, [inner_force - outer_force])
        norm_disp = __secant(x_value, fx_value)

    disp = __begetd(profile, norm_disp)
    return pressure, disp, inner_force, norm_disp
