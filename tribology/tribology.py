# -*- coding: utf-8 -*-

"""
Collection of tribology-related functions, developed and maintained at
KTH Royal Institute of Technology, Stockholm, Sweden.
"""

import copy
from math import log10 as lg
from math import sqrt, acos, cos, sin, pi, e

import numexpr as ne
import numpy as np
import scipy.sparse.linalg as spla


def secant_method(x_list, fx_list):
    """
    Applies secant method to find root x of function f(x).
    If not enough (x, f(x)) value pairs are known to apply
    secant method, a new x value is guessed by slightly changing the
    initial x value
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


def sliding_speed(vel_1, vel_2):
    """
    Calculate the sliding speed in a tribological contact based contact body
    velocities
    :param vel_1: velocity 1
    :param vel_2: velocity 2
    :return: sliding speed in contact between body 1 and 2
    """
    return vel_1 - vel_2


def rolling_speed(vel_1, vel_2):
    """
    Calculate the rolling speed in a tribological contact based contact body
    velocities
    :param vel_1:
    :param vel_2:
    :return:
    """
    return (vel_1 + vel_2) / 2


def srr(vel_1, vel_2):
    """
    Calculate the slide-to-roll ratio (srr) in a tribological contact
    :param vel_1: velocity 1
    :param vel_2: velocity 2
    :return: slide-to-roll ratio
    """
    return sliding_speed(vel_1, vel_2) / rolling_speed(vel_1, vel_2)


def rheometer_flat_torque(eta, gap_height, omega, r_a):
    """
    Calculate torque based on dynamic viscosity measurement from plate-on-
    plate rheometer
    :param eta: dynamic viscosity
    :param gap_height: distance between plates
    :param omega: angular velocity
    :param r_a: cone outer radius
    :return: torque
    """
    return eta * pi * omega * r_a ** 4 / (2 * gap_height)


def rheometer_flat_viscosity(torque, gap_height, omega, r_a):
    """
    Calculate the dynamic viscosity based on torque measurement from plate-on-
    plate rheometer
    :param torque: measured torque
    :param gap_height: distance between plates
    :param alpha: cone angle
    :param omega: angular velocity
    :param r_a: cone outer radius
    :return: dynamic viscosity
    """
    return 2 * gap_height * torque / (pi * omega * r_a ** 4)


def rheometer_cone_torque(eta, alpha, omega, r_a):
    """
    Calculate torque based on dynamic viscosity measurement from cone-on-
    plate rheometer
    :param eta: dinamic viscosity
    :param alpha: cone angle
    :param omega: angular velocity
    :param r_a: cone outer radius
    :return: torque
    """
    return 2 * pi * eta * omega * r_a ** 3 / (3 * alpha)


def rheometer_cone_viscosity(torque, alpha, omega, r_a):
    """
    Calculate the dynamic viscosity based on torque measurement from cone-on-
    plate rheometer
    :param torque: measured torque
    :param alpha: cone angle
    :param omega: angular velocity
    :param r_a: cone outer radius
    :return: dynamic viscosity
    """
    return 3 * torque * alpha / (2 * pi * omega * r_a ** 3)


def barus(eta_0, alpha_p, pressure):
    """
    Calculate dynamic viscosity at a given pressure based on the dynamic
    viscosity at ambient pressure and the pressure-viscosity coefficient
    :param eta_0: dynamic viscosity at ambient pressure
    :param alpha_p: pressure-viscosity coefficient
    :param pressure: pressure of interest
    :return: dynamic viscosity at pressure of interest
    """
    return eta_0 * e ** (alpha_p * pressure)


def rad_per_s_to_rpm(vel):
    """
    Convert from rpm to rad/s
    :param vel: velocity in rad/s
    :return: velocity in rpm
    """
    return 1 / rpm_to_rad_per_s(1 / vel)


def rpm_to_rad_per_s(vel):
    """
    Convert from rpm to rad/s
    :param vel: velocity in rpm
    :return: velocity in rad/s
    """
    return vel / 60 * 2 * pi


def ball_on_3_plates_slide_radius(r_ball, plate_angle=1.5708):
    """
    Sliding radius (lever arm) for ball-on-3-plates setup
    :param r_ball: radius of rotating ball
    :param plate_angle: plate angle with respect to ball in rad
                        (default corresponds to 45 degree)
    :return: sliding radius
    """
    return r_ball * sin((pi - plate_angle) / 2)


def ball_on_3_plates_norm_force(ax_force, plate_angle=1.5708):
    """
    Calculate normal force per contact in ball-on-3-plates setup
    :param ax_force: axial force on rotating ball
    :param plate_angle: plate angle with respect to each other
                        (default corresponds to 90 degree)
    :return: normal force per contact
    """
    return ax_force / 3 / cos(plate_angle / 2)


def fourball_geometry(r_1, r_2):
    """
    Geometric parameters of 4-ball setup
    :param r_1: radius rotating ball
    :param r_2: radius stationary balls
    :return: sliding radius (lever arm) on rotating ball, contact angle in rad
    """
    r_circum_circle = sqrt(3) / 3 * 2 * r_2
    contact_angle = acos(r_circum_circle / (r_1 + r_2))
    sliding_radius = r_circum_circle - r_2 * cos(contact_angle)
    return sliding_radius, contact_angle


def fourball_norm_force(r_1, r_2, ax_force):
    """
    Calculate normal force per contact in 4 ball setup
    :param r_1: radius rotating ball
    :param r_2: radius stationary balls
    :param ax_force: axial force on rotating ball
    :return: normal force per contact
    """
    _, contact_angle = fourball_geometry(r_1, r_2)
    return ax_force / sin(contact_angle) / 3


def convert_prefix(value, in_prefix, out_prefix):
    """
    Convert between different SI unit prefixes
    :param value:  value to convert
    :param in_prefix:  input SI unit prefix
    :param out_prefix:  output SI unit prefix
    :return:  value converted to new unit prefix
    """
    prefix = {'p': 10 ** -12,
              'n': 10 ** -9,
              'mu': 10 ** -6,
              'm': 10 ** -3,
              '': 10 ** -0,
              'k': 10 ** 3,
              'M': 10 ** 6,
              'G': 10 ** 9,
              'T': 10 ** 12}
    return value * prefix[in_prefix] / prefix[out_prefix]


def abbott_firestone(trace, res=100):
    """
    Calculate Abbott-Firestone curve for 2D surface trace
    :param trace: vector containing surface heights
    :param res: number of discontinuous data steps
    :return: data baskets (steps) and cumulative distribution
             (x and y data for Abbott-Firestone plot)
    """
    baskets = np.linspace(np.amax(trace), np.amin(trace), res)
    cum_dist = []
    for basket in baskets:
        cum_dist.append(len(np.where(trace >= basket)[0]) / len(trace) * res)
    return baskets, cum_dist


def kin2dyn(kin, density):
    """
    kinematic to dynamic viscosity
    :param kin: kinematic viscosity
    :param density: density
    :return: dynamic viscosity
    """
    return kin / density


def dyn2kin(dyn, density):
    """
    dynamic to kinematic viscosity
    :param dyn: dynamic viscosity
    :param density: density
    :return: kinematic viscosity
    """
    return dyn * density





def walther_calc_z_from_nu(kin_visc):
    """
    calculate z factor for Walther Equation (ASTM D341)
    :param kin_visc: kinematic viscosity
    :return: z factor
    """
    return kin_visc + 0.7 + \
           10 ** (-1.47 - 1.84 * kin_visc - 0.51 * kin_visc ** 2)


def walther_calc_nu(zed):
    """
    calculate kinematic viscosity for Walther equation (ASTM D341)
    :param zed: z factor
    :return: kinematic viscosity
    """
    return (zed - 0.7) - 10 ** (-0.7487 - 3.295 * (zed - 0.7) +
                                0.6119 * (zed - 0.7) ** 2 - 0.3193 *
                                (zed - 0.7) ** 3)


def walther(temp_1, nu_1, temp_2, nu_2, temp_3):
    """
    Calculates kinematic viscosity at temperature temp_3 based on kin.
    viscosities at temperatures temp_1 and tempo_2. Equations according to
    ASTM D341.
    :param temp_1: temperature in deg C
    :param nu_1: kin. viscosity at t1 in cSt
    :param temp_2: temperature in deg C
    :param nu_2: kin. viscosity at t2 in cSt
    :param temp_3: temperature of interest, in deg C
    :return: kin. viscosity in cSt at temperature t3
    """
    abs_zero = -273.15
    viscs = [nu_1, nu_2]
    thetas = [temp_1 - abs_zero, temp_2 - abs_zero, temp_3 - abs_zero]

    zed = [walther_calc_z_from_nu(nu) for nu in viscs]
    const_a = (lg(lg(zed[0])) - lg(lg(zed[1])) *
               lg(thetas[0]) / lg(thetas[1])) / \
              (1 - lg(thetas[0]) / lg(thetas[1]))
    const_b = (const_a - lg(lg(zed[1]))) / lg(thetas[1])
    return walther_calc_nu(10 ** 10 ** (const_a - const_b * lg(thetas[2])))


def effective_radius(r_1, r_2):
    """
    Effective radius according to Hertzian contact theory
    :param r_1: radius 1
    :param r_2: radius 2
    :return: reduced (effective) radius
    """
    if (r_1 == 0 or abs(r_1) == float('inf')) and r_2 != 0:
        return r_2
    elif (r_2 == 0 or abs(r_2) == float('inf')) and r_1 != 0:
        return r_1
    elif (r_1 == 0 or abs(r_1) == float('inf')) and \
            (r_2 == 0 or abs(r_2) == float('inf')):
        return 0
    elif r_1 == -r_2:
        return 0
    else:
        return 1 / (1 / r_1 + 1 / r_2)


def effective_radii(r_x_1, r_y_1, r_x_2, r_y_2, ):
    """
    Effective radii for combination of 2 bodies according to Hertzian
    contact theory
    :param r_x_1: radius of body 1 in x
    :param r_y_1: radius of body 1 in y
    :param r_x_2: radius of body 2 in x
    :param r_y_2: radius of body 2 in y
    :return: effective radius (total) and effective radii in each plane
    """
    recip_radius = []
    for radius in [r_x_1, r_y_1, r_x_2, r_y_2]:
        if radius == 0:
            recip_radius.append(float('inf'))
        else:
            recip_radius.append(1 / radius)

    r_eff_x = effective_radius(r_x_1, r_x_2)
    r_eff_y = effective_radius(r_y_1, r_y_2)
    return effective_radius(r_eff_x, r_eff_y), r_eff_x, r_eff_y


def effective_modulus(e_1, nu_1, e_2, nu_2):
    """
    Effective Young's modulus according to Hertzian contact theory
    :param e_1: young's modulus body 1
    :param nu_1: poisson ratio body 1
    :param e_2: young's modulus body 2
    :param nu_2: poisson ratio body 2
    :return:
    """
    return 1 / ((1 - nu_1 ** 2) / (2 * e_1) + (1 - nu_2 ** 2) / (2 * e_2))


def influ_mat_coordinate_grid(axis):
    """
    Generate coordinate grid based on axis as required for influence matrix
    generation
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


def influ_mat_reduce(influ_mat):
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


def boundary_element_influ_mat(x_axis, y_axis, e_eff):
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
    a_factor = (x_axis[-1] - x_axis[0]) / len_x / 2
    b_factor = y_axis[-1] - y_axis[0] / len_x / 2
    x_grid = influ_mat_coordinate_grid(x_axis)
    y_grid = influ_mat_coordinate_grid(y_axis)

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

    return influence_matrix_complete * e_eff / pi


def get_displacements(profile, norm_disp):
    """
    Calculate local elastic displacements of profile
    :param profile: combined profile of two contacting bodies
    :param norm_disp: global normal elastic displ between contacting bodies
    :return: local normal displacements as a result of global displ norm_displ
    """
    displ_field = np.subtract(np.ones(profile.shape) * norm_disp, profile)
    disp = np.reshape(displ_field, (profile.shape[0] * profile.shape[1]))
    return disp


def combine_profile(profile_1, profile_2):
    """
    Combine two body profiles for boundary element calculation
    :param profile_1: array-like profile heights
    :param profile_2: array-like profile heights
    :return: negative combined profile heights
    """
    return -(profile_1 + profile_2)


def boundary_element_solve_pressure(profile_1, profile_2, outer_force,
                                    red_influ_mat, delta_x, delta_y,
                                    norm_disp=0.001, max_offset=0.005):
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
    profile = combine_profile(profile_1, profile_2)

    # while difference between inner forces and outer forces is significant
    while abs(inner_force - outer_force) > max_offset * outer_force:

        # update local displacements
        disp = get_displacements(profile, norm_disp)

        # find negative pressure arguments
        pressure = spla.gmres(red_influ_mat, disp)[0]
        p_index = np.zeros(len(pressure))
        negative_p = np.where(pressure < 0)[0]
        p_neg = copy.deepcopy(negative_p)

        # remove elements with negative pressure
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

        # calculate resulting force and adjust displacement for next loop
        pressure = np.reshape(pressure, (profile.shape[0], profile.shape[1]))
        inner_force = sum(sum(np.multiply(delta_x * delta_y, pressure)))
        x_value = np.append(x_value, [norm_disp])
        fx_value = np.append(fx_value, [inner_force - outer_force])
        norm_disp = secant_method(x_value, fx_value)

    disp = get_displacements(profile, norm_disp)
    return pressure, disp, inner_force, x_value[-1]



if __name__ == "__main__":
    pass
