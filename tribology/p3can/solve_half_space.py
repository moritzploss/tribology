import copy
import math
import os
import uuid

import numexpr as ne
import numpy as np
import scipy.sparse.linalg as spla

from Constants import PrintOpts, SubDir
from numeric_methods import secant_method
from system_functions import print_progress, print_it, make_directory, \
    print_rel_prox_to_sol


def solve_for_pressure(profile1, profile2, normal_load, red_influ_mat,
                       init_displ, delta_x, delta_y, print_prog=True):
    """
    Solve system of equations:
    
        [pressure] = [influence matrix]^-1 * [displacement]
        
    Stop solver once inner and outer forces are in equilibrium
    """

    # initialise variables
    x_value = [0]
    fx_value = [normal_load]
    pressure = []
    resulting_force = 0
    profile = -(profile1 + profile2)
    threshold_factor = 0.005

    # do while difference between inner forces and outer forces is significant
    while abs(resulting_force - normal_load) > threshold_factor * normal_load:
        displacement_field = np.subtract(np.ones(profile.shape) * init_displ,
                                         profile)
        u = np.reshape(displacement_field,
                       (profile.shape[0] * profile.shape[1]))

        # lots of options to solve this, the active implementation proved to be
        # the fastest for me
        # p = np.linalg.solve(red_influ_mat, u)
        p = spla.gmres(red_influ_mat, u)[0]

        p_index = np.zeros(len(p))
        negative_p = np.where(p < 0)[0]
        p_neg = copy.deepcopy(negative_p)

        # do while there are elements with negative pressure
        while len(negative_p) > 0:
            p[p_neg] = 0
            p_index[p_neg] = 1
            u_new_reduced = np.delete(u, [p_neg], axis=0)
            g_new_reduced = np.delete(red_influ_mat, [p_neg], axis=0)
            g_new_reduced = np.delete(g_new_reduced, [p_neg], axis=1)
            if p[np.where(p_index == 0)].size > 0:
                # again, lots of options available here, this seems to be the
                # fastest for most cases
                p[np.where(p_index == 0)] = \
                spla.gmres(g_new_reduced, u_new_reduced)[0]
            negative_p = np.where(p < 0)[0]
            p_neg = np.append(p_neg, negative_p)

        # calculate resulting force and adjust displacement for next calculation
        # loop
        pressure = np.reshape(p, (profile.shape[0], profile.shape[1]))
        resulting_force = sum(sum(np.multiply(delta_x * delta_y, pressure)))
        x_value = np.append(x_value, [init_displ])
        fx_value = np.append(fx_value, [resulting_force - normal_load])
        init_displ = secant_method(x_value, fx_value)
        if print_prog is True:
            print_rel_prox_to_sol(resulting_force, normal_load,
                                  threshold_factor)
    return pressure, resulting_force, x_value[-1]


def pre_solve_half_space(profile1, profile2, x_axis, y_axis, res_x, res_y,
                         delta_x, delta_y, e1, e2, ny1, ny2,
                         normal_load, init_displ=0.0005, print_prog=False):
    """Solve half space with simplified geometries to find grid size for actual
    simulation"""
    red_influ_mat = calc_reduced_influence_matrix(x_axis, y_axis, res_x, res_y,
                                                  delta_x, delta_y, e1, e2, ny1,
                                                  ny2,
                                                  print_prog=False)
    pressure, resulting_force, init_displ = solve_for_pressure(profile1,
                                                               profile2,
                                                               normal_load,
                                                               red_influ_mat,
                                                               init_displ,
                                                               delta_x, delta_y,
                                                               print_prog)
    return pressure, init_displ


def solve_half_space(profile1, profile2, x_axis, y_axis, res_x, res_y, delta_x,
                     delta_y, e1, e2, ny1, ny2, normal_load,
                     res_dir, influ_mat_db='empty', print_prog=True,
                     init_displ=0.0005):
    """Orchestrate half space calculation and influence matrix
    generation/retrieval"""
    influ_mat_path = os.sep.join([res_dir, SubDir.np_db.value, influ_mat_db])
    if os.path.exists(influ_mat_path):
        red_influ_mat = np.load(influ_mat_path)
    else:
        print_it('calculating influence matrix', PrintOpts.lvl1.value)
        red_influ_mat = calc_reduced_influence_matrix(x_axis, y_axis, res_x,
                                                      res_y, delta_x, delta_y,
                                                      e1, e2, ny1, ny2,
                                                      print_prog=True)
        print_it('saving influence matrix', PrintOpts.lvl1.value)
        np_db_directory = make_directory(res_dir, SubDir.np_db.value)
        influ_mat_db = 'influence-matrix-{}.npy'.format(uuid.uuid4())
        np.save(os.sep.join([np_db_directory, influ_mat_db]), red_influ_mat)
        print_it('solving first half space', PrintOpts.lvl1.value)

    pressure, resulting_force, init_displ = \
        solve_for_pressure(profile1, profile2, normal_load, red_influ_mat,
                           init_displ, delta_x, delta_y, print_prog)
    return pressure, influ_mat_db


def calc_reduced_influence_matrix(x_axis, y_axis, res_x, res_y, delta_x,
                                  delta_y, e1, e2, ny1, ny2, print_prog=False):
    """Calculate the complete influence matrix C and extract a reduced influence
    matrix"""

    # initialise variables
    influence_matrix_complete = np.zeros((res_x, res_y, res_x, res_y))
    reduced_influence_matrix = np.zeros((res_x * res_y, res_x * res_y))
    youngs_mod_over_pi = ((1 - math.pow(ny1, 2)) / e1 + (
        1 - math.pow(ny2, 2)) / e2) / math.pi
    a = delta_x / 2
    b = delta_y / 2

    # generate x coordinate grid
    x_new = np.zeros((1, res_x))
    x_new[0, :] = x_axis
    vertical_x_axis = np.zeros((res_x, 1))
    vertical_x_axis[:, 0] = x_axis
    x = np.repeat(x_new, res_x, axis=0)
    x = np.absolute(np.subtract(x, vertical_x_axis))

    # generate y coordinate grid
    y_new = np.zeros((1, res_y))
    y_new[0, :] = y_axis
    vertical_y_axis = np.zeros((res_y, 1))
    vertical_y_axis[:, 0] = y_axis
    y = np.repeat(y_new, res_y, axis=0)
    y = np.absolute(np.subtract(y, vertical_y_axis))

    # use numexpr to evaluate expressions
    xpa = ne.evaluate('x + a')
    xma = ne.evaluate('x - a')
    ypb = ne.evaluate('y + b')
    ymb = ne.evaluate('y - b')

    '''
    The code in comments below is a simple implementation to calculate the 
    reduced influence matrix. Due to its 4 nested for-loops, it's quite slow, 
    but it's readable. The active implementation (see further below) is 
    basically the same code, but in vectorized from, which makes it much faster, 
    but less readable. It is possible to vectorize the active implementation 
    even further and to get rid of another two for-loops; however, this leads 
    to matrix operations on very large matrices, and my notebook had problems 
    with that at larger grid sizes. also, for high grid resolution
    the current bottleneck is solving the half space, not generating the 
    influence matrix.
        -- Moritz
    
    influence_matrix = np.zeros((res_x, res_y))
    counter = 0
    for i in range(0, res_x):
        print_progress(i, res_x)
        for j in range(0, res_y):
            for i_prime in range(0, res_x):
                x = abs(x_axis[i] - x_axis[i_prime])
                for j_prime in range(0, res_y):
                    y = abs(y_axis[j] - y_axis[j_prime])

                    influence_matrix[i_prime, j_prime] = youngs_mod_over_pi * \
                            (
                        (x + a) * math.log(((y + b) + math.sqrt(math.pow((y + b), 2) + math.pow((x + a), 2))) / (
                            (y - b) + math.sqrt(math.pow((y - b), 2) + math.pow((x + a), 2)))) +
                        (y + b) * math.log(((x + a) + math.sqrt(math.pow((y + b), 2) + math.pow((x + a), 2))) / (
                            (x - a) + math.sqrt(math.pow((y + b), 2) + math.pow((x - a), 2)))) +
                        (x - a) * math.log(((y - b) + math.sqrt(math.pow((y - b), 2) + math.pow((x - a), 2))) / (
                            (y + b) + math.sqrt(math.pow((y + b), 2) + math.pow((x - a), 2)))) +
                        (y - b) * math.log(((x - a) + math.sqrt(math.pow((y - b), 2) + math.pow((x - a), 2))) / (
                            (x + a) + math.sqrt(math.pow((y - b), 2) + math.pow((x + a), 2))))
                            )

                    #influence_matrix[i_prime, j_prime] = ne.evaluate(eval_string)
            reduced_influence_matrix[counter, :] = np.reshape(influence_matrix, (res_x * res_y))
            counter += 1
    '''

    # this code is a vectorized version of the above code. the code assembles
    # the complete influence matrix
    for j in range(0, res_y):
        if print_prog is True:
            print_progress(j, res_y)
        for j_prime in range(0, res_y):
            influence_matrix_complete[:, j, :, j_prime] = \
                (
                    np.multiply(xpa, np.log(
                        np.divide(
                            ((ypb[j, j_prime]) + np.sqrt(
                                np.multiply((ypb[j, j_prime]),
                                            (ypb[j, j_prime])) +
                                np.multiply(xpa, xpa)))
                            ,
                            ((ymb[j, j_prime]) + np.sqrt(
                                np.multiply((ymb[j, j_prime]),
                                            (ymb[j, j_prime])) +
                                np.multiply(xpa, xpa)))
                        )
                    ))
                    +
                    (ypb[j, j_prime]) * np.log(
                        np.divide(
                            (xpa + np.sqrt(np.multiply((ypb[j, j_prime]),
                                                       (ypb[j, j_prime])) +
                                           np.multiply(xpa, xpa)))
                            ,
                            (xma + np.sqrt(np.multiply((ypb[j, j_prime]),
                                                       (ypb[j, j_prime])) +
                                           np.multiply(xma, xma)))
                        )
                    )
                    +
                    np.multiply(xma, np.log(
                        np.divide(
                            ((ymb[j, j_prime]) + np.sqrt(
                                np.multiply((ymb[j, j_prime]),
                                            (ymb[j, j_prime])) +
                                np.multiply(xma, xma)))
                            ,
                            ((ypb[j, j_prime]) + np.sqrt(
                                np.multiply((ypb[j, j_prime]),
                                            (ypb[j, j_prime])) +
                                np.multiply(xma, xma)))
                        )
                    ))
                    +
                    (ymb[j, j_prime]) * np.log(
                        np.divide(
                            (xma + np.sqrt(np.multiply((ymb[j, j_prime]),
                                                       (ymb[j, j_prime])) +
                                           np.multiply(xma, xma)))
                            ,
                            (xpa + np.sqrt(np.multiply((ymb[j, j_prime]),
                                                       (ymb[j, j_prime])) +
                                           np.multiply(xpa, xpa)))
                        )
                    )
                )

    influence_matrix_complete *= youngs_mod_over_pi

    # extract the reduced influence matrix from the complete influence matrix
    counter = 0
    for i in range(0, res_x):
        for j in range(0, res_y):
            reduced_influence_matrix[counter, :] = np.reshape(
                influence_matrix_complete[i, j, :, :], (res_x * res_y))
            counter += 1

    return reduced_influence_matrix
