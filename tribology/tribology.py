# -*- coding: utf-8 -*-

"""
Methods related to general tribology
"""

from math import sqrt, acos, cos, sin, pi, floor

import numpy as np


def profball(x_axis, r_ball):
    """
    Generate a 2D ball profile based on ball radius and ball axis
    :param x_axis: vector containing coordinate points for which to calculate
                   profile heights
    :param r_ball: the radius of the ball
    :return:
    """
    prof = np.array(abs(r_ball) - np.sqrt(r_ball ** 2 - np.power(x_axis, 2)))
    return prof * np.sign(r_ball)


def profrevolve(prof_2d, y_axis, y_diam):
    """
    creates a 3d profile by revolving a 2D profile around the central axis of
    a body
    :param prof_2d: 2d profile vector containing profile heights
    :param y_axis: vector containing coordinate points for which to calculate
                   profile heights in y-direction
    :param y_diam: diameter around which to revolve prof_2d
    :return:
    """
    len_x = len(prof_2d)
    len_y = len(y_axis)
    prof_3d = np.zeros((len_x, len_y))
    sign_diam = np.sign(y_diam)
    for i_x in range(len_x):
        for i_y in range(len_y):
            bracket = pow((y_diam / 2 - sign_diam * prof_2d[i_x]), 2)\
                      - pow(y_axis[i_y], 2)
            if bracket <= 0:
                prof_3d[i_x, i_y] = abs(y_diam) / 2
            else:
                prof_3d[i_x, i_y] = abs(y_diam) / 2 - sign_diam * sqrt(bracket)
    if y_diam > 0:
        min_prof = prof_3d.min()
    else:
        min_prof = prof_3d[round(len_x / 2) - 1, round(len_x / 2) - 1]
    prof_3d = -(prof_3d - min_prof)
    y_profile = prof_3d[:, floor(len_y / 2)]
    return -prof_3d, y_profile


def vslide(vel_1, vel_2):
    """
    Calculate the sliding speed in a tribological contact based contact body
    velocities
    :param vel_1: velocity 1
    :param vel_2: velocity 2
    :return: sliding speed in contact between body 1 and 2
    """
    return vel_1 - vel_2


def vroll(vel_1, vel_2):
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
    return vslide(vel_1, vel_2) / vroll(vel_1, vel_2)


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


def rball3plates(r_ball, plate_angle=1.5708):
    """
    Sliding radius (lever arm) for ball-on-3-plates setup
    :param r_ball: radius of rotating ball
    :param plate_angle: plate angle with respect to ball in rad
                        (default corresponds to 45 degree)
    :return: sliding radius
    """
    return r_ball * sin((pi - plate_angle) / 2)


def fball3plates(ax_force, plate_angle=1.5708):
    """
    Calculate normal force per contact in ball-on-3-plates setup
    :param ax_force: axial force on rotating ball
    :param plate_angle: plate angle with respect to each other
                        (default corresponds to 90 degree)
    :return: normal force per contact
    """
    return ax_force / 3 / cos(plate_angle / 2)


def gfourball(r_1, r_2):
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


def ffourball(r_1, r_2, ax_force):
    """
    Calculate normal force per contact in 4 ball setup
    :param r_1: radius rotating ball
    :param r_2: radius stationary balls
    :param ax_force: axial force on rotating ball
    :return: normal force per contact
    """
    _, contact_angle = gfourball(r_1, r_2)
    return ax_force / sin(contact_angle) / 3


def convpref(value, p_in='', p_out=''):
    """
    Convert between different SI unit prefixes
    :param value:  value to convert
    :param p_in:  input SI unit prefix
    :param p_out:  output SI unit prefix
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
    return value * prefix[p_in] / prefix[p_out]


def abbottfirestone(trace, res=100):
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


def reff(r_1, r_2):
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


def eeff(r_x_1, r_y_1, r_x_2, r_y_2, ):
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

    r_eff_x = reff(r_x_1, r_x_2)
    r_eff_y = reff(r_y_1, r_y_2)
    return reff(r_eff_x, r_eff_y), r_eff_x, r_eff_y


def meff(e_1, nu_1, e_2, nu_2):
    """
    Effective Young's modulus according to Hertzian contact theory
    :param e_1: young's modulus body 1
    :param nu_1: poisson ratio body 1
    :param e_2: young's modulus body 2
    :param nu_2: poisson ratio body 2
    :return:
    """
    return 1 / ((1 - nu_1 ** 2) / (2 * e_1) + (1 - nu_2 ** 2) / (2 * e_2))


if __name__ == "__main__":
    pass
