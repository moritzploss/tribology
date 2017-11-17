# -*- coding: utf-8 -*-

"""

Methods related to boundary element solvers for contact mechanics calculations

"""


import copy
from math import pi

import numexpr as ne
import numpy as np
import scipy.sparse.linalg as spla


def secant(x_list, fx_list):
    """

    Applies secant method to find root *x* of function *f(x)*.
    If not enough (*x*, *f(x)*) value pairs are available to apply
    secant method, a new *x* value is guessed by slightly changing the
    initial *x* value

    :param x_list: list of function inputs
    :param fx_list: list of function outputs

    :return: estimate of function root

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

    Private helper method for influence matrix generation

    :param axis: axis to make grid

    :return: grid

    """
    len_axis = len(axis)
    vec = np.zeros((1, len_axis))
    vec[0, :] = axis
    vertical_ax = np.zeros((len_axis, 1))
    vertical_ax[:, 0] = axis
    grid = np.repeat(vec, len_axis, axis=0)
    return np.absolute(np.subtract(grid, vertical_ax))


def beinflumatred(influ_mat):
    """

    Extract the reduced influence matrix from the complete influence matrix

    :param influ_mat: complete influence matrix

    :return: reduced influence matrix

    """
    shape_mat = np.shape(influ_mat)
    len_mat = shape_mat[0] * shape_mat[1]
    reduced_influence_matrix = np.zeros((len_mat, len_mat))
    counter = 0
    for i in range(0, shape_mat[0]):
        for j in range(0, shape_mat[1]):
            reduced_influence_matrix[counter, :] = \
                np.reshape(influ_mat[i, j, :, :], len_mat)
            counter += 1
    return reduced_influence_matrix


def beinflumat(x_axis, y_axis, e_eff):
    """

    Generate an influence matrix as required for boundary element contact
    mechanics calculations

    :param x_axis: x-axis of coordinate grid
    :param y_axis: y-axis of coordinate grid
    :param e_eff: effective youngs modulus

    :return: complete influence matrix

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

    return influence_matrix_complete * 1 / e_eff * 2 / pi


def __begetd(profile, norm_disp):
    """

    Calculate local elastic displacements of profile

    :param profile: combined profile of two contacting bodies
    :param norm_disp: global normal elastic displ between contacting bodies

    :return: local normal displacements as a result of global displ norm_displ

    """
    displ_field = np.subtract(np.ones(profile.shape) * norm_disp, profile)
    disp = np.reshape(displ_field, (profile.shape[0] * profile.shape[1]))
    return disp


def combineprof(profile_1, profile_2):
    """

    Combine two body profiles for boundary element calculation

    :param profile_1: array-like profile heights
    :param profile_2: array-like profile heights

    :return: negative combined profile heights

    """
    return profile_1 + profile_2


def __solvepress(red_influ_mat, disp):
    """

    solve for pressure distribution

    :param red_influ_mat: reduced influence matrix
    :param disp: displacement field

    :return: pressure field

    """

    # find negative pressure arguments
    pressure = spla.gmres(red_influ_mat, disp)[0]
    p_index = np.zeros(len(pressure))
    negative_p = np.where(pressure < 0)[0]
    p_neg = copy.deepcopy(negative_p)

    while len(negative_p) > 0:
        pressure[p_neg] = 0
        p_index[p_neg] = 1
        u_new_reduced = np.delete(disp, [p_neg], axis=0)
        g_new_reduced = np.delete(red_influ_mat, [p_neg], axis=0)
        g_new_reduced = np.delete(g_new_reduced, [p_neg], axis=1)
        if pressure[np.where(p_index == 0)].size > 0:
            pressure[np.where(p_index == 0)] = \
                spla.gmres(g_new_reduced, u_new_reduced)[0]
        negative_p = np.where(pressure < 0)[0]
        p_neg = np.append(p_neg, negative_p)

    return pressure


def besolve(profile_1, profile_2, outer_force,
            red_influ_mat, delta_x, delta_y,
            norm_disp=0.1, max_offset=0.005):
    """

    Solve system of equations:

        [pressure] = [influence matrix]^-1 * [displacement]

    Stop solver once inner and outer forces are in equilibrium

    :param profile_1: array with profile heights of body 1
    :param profile_2: array with profile hights of body 2
    :param outer_force: outer (normal) force on body 1/2
    :param red_influ_mat: reduced influence matrix
    :param delta_x: grid spacing of profile arrays in x-direction
    :param delta_y: grid spacing of profile arrays in y-direction
    :param norm_disp: initial normal elastic deformation to start calculation
    :param max_offset: maximum allowed percentage difference of inner and outer
                       force at end of calculation

    :return: local pressure (array), local displacements (array).
             inner force (scalar), global normal displacement (scalar)

    """
    # initialise variables
    x_value = [0]
    fx_value = [outer_force]
    pressure = 0
    inner_force = 0
    profile = combineprof(profile_1, profile_2)

    # while difference between inner forces and outer forces is significant
    while abs(inner_force - outer_force) > max_offset * outer_force:

        # update local displacements
        disp = __begetd(profile, norm_disp)
        pressure = __solvepress(red_influ_mat, disp)

        # calculate resulting force and adjust displacement for next loop
        pressure = np.reshape(pressure, (profile.shape[0], profile.shape[1]))
        inner_force = sum(sum(np.multiply(delta_x * delta_y, pressure)))
        x_value = np.append(x_value, [norm_disp])
        fx_value = np.append(fx_value, [inner_force - outer_force])
        norm_disp = secant(x_value, fx_value)

    disp = __begetd(profile, norm_disp)
    return pressure, disp, inner_force, x_value[-1]
