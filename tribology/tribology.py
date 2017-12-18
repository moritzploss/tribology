# -*- coding: utf-8 -*-

"""

The functions in this module are related to tribology in general.

"""

from math import sqrt, acos, cos, sin, pi, floor

import numpy as np


def profrolleriso(x_axis, diam, length):
    """

    Generate a roller profile according to DIN 26281 for rollers with
    :code:`length` :math:`\\leq` :code:`2.5` :math:`\\cdot` :code:`diam`.

    Parameters
    ----------
    x_axis: ndarray
        The x-axis of the roller. The origin of the axis must be at the center
        of the roller. To obtain sensible results, the axis limits must be
        within :code:`+/- length / 2`.
    diam: scalar
        The diameter of the roller.
    length: scalar
        The length of the roller.

    Returns
    -------
    x_profile: ndarray
        The profile heights of the roller along :code:`x_axis`.

    """
    len_ax = abs(max(x_axis) - min(x_axis))
    if len_ax == length:
        ax_x = x_axis[1:-1]
    else:
        ax_x = x_axis

    prof = 0.00035 * diam * np.log(1 / (1 - np.power(((2 * ax_x) / length), 2)))
    prof += np.min(prof)

    if len_ax == length:
        diam_arr = np.asarray([diam])
        x_profile = np.concatenate((diam_arr, np.concatenate((prof, diam_arr))))
    else:
        x_profile = prof

    return x_profile


def profball(x_axis, r_ball):
    """

    Create a 1D ball profile along a central axis.

    Parameters
    ----------
    x_axis: ndarray
        The coordinate points for which to calculate profile heights. Minimum
        and maximum value must be within +/- **r_ball**.
    r_ball: scalar
        The radius of the ball.

    Returns
    -------
    prof_heights: ndarray
        The profile heights along the **x_axis**.

    """
    prof = np.array(abs(r_ball) - np.sqrt(r_ball ** 2 - np.power(x_axis, 2)))
    prof_heights = prof * np.sign(r_ball)
    return prof_heights


def profrevolve(prof_1d, y_axis, y_diam):
    """

    Create a 2D surface (profile heights) by revolving a 1D profile around a
    central axis.

    Parameters
    ----------
    prof_1d: ndarray
        The profile heights along the central (x) axis of the surface.
    y_axis: ndarray
        The y-axis of the surface.
    y_diam: scalar
        The maximum diameter of the surface along its y-axis.

    Returns
    -------
    prof_2d: ndarray
        The surface heights of the profile. The array size is
        :code:`len(prof_2d)` :math:`\\times` :code:`len(y_axis)`
    y_profile: ndarray
        The profile heights along the y-direction of the surface at the central
        element of the x-axis.

    """
    len_x = len(prof_1d)
    len_y = len(y_axis)
    prof_2d = np.zeros((len_x, len_y))
    sign_diam = np.sign(y_diam)
    for i_x in range(len_x):
        for i_y in range(len_y):
            bracket = pow((y_diam / 2 - sign_diam * prof_1d[i_x]), 2)\
                      - pow(y_axis[i_y], 2)
            if bracket <= 0:
                prof_2d[i_x, i_y] = abs(y_diam) / 2
            else:
                prof_2d[i_x, i_y] = abs(y_diam) / 2 - sign_diam * sqrt(bracket)
    if y_diam > 0:
        min_prof = prof_2d.min()
    else:
        min_prof = prof_2d[round(len_x / 2) - 1, round(len_x / 2) - 1]
    prof_2d = (prof_2d - min_prof)
    y_profile = prof_2d[:, floor(len_y / 2)]
    return prof_2d, y_profile


def vslide(vel_1, vel_2):
    """

    Calculate the sliding speed in a tribological contact based on contact body
    velocities.

    Parameters
    ----------
    vel_1: ndarray, scalar
        The contact velocity of body 1.
    vel_2: ndarray, scalar
        The contact velocity of body 2

    Returns
    -------
    vel_slide: ndarray, scalar
        The sliding velocity in the tribological contact.

    """
    vel_slide = vel_1 - vel_2
    return vel_slide


def vroll(vel_1, vel_2):
    """

    Calculate the rolling speed in a tribological contact based on contact body
    velocities.

    Parameters
    ----------
    vel_1: ndarray, scalar
        The contact velocity of body 1.
    vel_2: ndarray, scalar
        The contact velocity of body 2

    Returns
    -------
    vel_roll: ndarray, scalar
        The rolling velocity in the tribological contact.

    """
    vel_roll = (vel_1 + vel_2) / 2
    return vel_roll


def srr(vel_1, vel_2):
    """

    Calculate the slide-to-roll ratio (SRR) in a tribological contact based on
    contact body velocities.

    Parameters
    ----------
    vel_1: ndarray, scalar
        The contact velocity of body 1.
    vel_2: ndarray, scalar
        The contact velocity of body 2

    Returns
    -------
    srr: ndarray, scalar
        The slide-to-roll ratio in the tribological contact.

    """
    val_srr = vslide(vel_1, vel_2) / vroll(vel_1, vel_2)
    return val_srr


def radpersec2rpm(vel_rad_per_sec):
    """

    Convert velocity from rotations per minute (rpm) to radians per second
    (rad/s).


    Parameters
    ----------
    vel_rad_per_sec: ndarray, scalar
        The velcity in radians per second.

    Returns
    -------
    vel_rpm: ndarray, scalar
        The velocity in rotations per minute.

    """
    vel_rpm = 1 / rpm2radpersec(1 / vel_rad_per_sec)
    return vel_rpm


def rpm2radpersec(vel_rpm):
    """

     Convert velocity from radians per second (rad/s) to rotations per minute
     (rpm).


     Parameters
     ----------
     vel_rpm: ndarray, scalar
         The velcity in rotations per minute.

     Returns
     -------
     vel_rad_per_sec: ndarray, scalar
         The velocity in radians per second.

     """
    vel_rad_per_sec = vel_rpm / 60 * 2 * pi
    return vel_rad_per_sec


def rball3plates(r_ball, plate_angle=1.5708):
    """

    Calculate the sliding radius (lever arm) for a ball-on-3-plates test setup.

    Parameters
    ----------
    r_ball: ndarray, scalar
        The radius of the ball.
    plate_angle: ndarray, scalar, optional
        The plate angle in radians, measured with respect to the rotational
        axis of the ball. Default value corresponds to 45 degree.

    Returns
    -------
    r_slide: ndarray, scalar
        The sliding radius.

    """
    r_slide = r_ball * np.sin((pi - plate_angle) / 2)
    return r_slide


def fball3plates(ax_force, plate_angle=1.5708):
    """

    Calculate the normal force per contact in ball-on-3-plates setup.

    Parameters
    ----------
    ax_force: ndarray, scalar
        The force acting along the rotational axis of the ball (axial force).
    plate_angle: ndarray, scalar, optional
        The plate angle in radians, measured with respect to the rotational
        axis of the ball. Default value corresponds to 45 degree.

    Returns
    -------
    norm_force:
        The normal force acting in each ball-plate contact.

    """
    norm_force = ax_force / 3 / np.cos(plate_angle / 2)
    return norm_force


def gfourball(r_1, r_2):
    """

    Geometric parameters of 4-ball setup

    Parameters
    ----------
    r_1 : scalar
        The radius of the rotating ball.
    r_2 : scalar
        The radius of the (a) stationary ball.

    Returns
    -------
    sliding_radius : scalar
        The sliding radius (lever arm) on the rotating ball.
    contact_angle : scalar
        Contact angle between direction of ball normal force and vertical axis.

    """
    r_circum_circle = sqrt(3) / 3 * 2 * r_2
    contact_angle = acos(r_circum_circle / (r_1 + r_2))
    sliding_radius = r_circum_circle - r_2 * cos(contact_angle)
    return sliding_radius, contact_angle


def ffourball(r_1, r_2, ax_force):
    """

    Calculate the normal force per contact in a 4-ball test setup.

    Parameters
    ----------
    r_1:
        The radius of the rotating ball.
    r_2:
        The radius of the stationary balls.
    ax_force:
        The force acting along the rotational axis of the rotating ball.

    Returns
    -------
    norm_force: ndarray, scalar
        The normal force in a single ball-ball contact.

    """
    _, contact_angle = gfourball(r_1, r_2)
    norm_force = ax_force / sin(contact_angle) / 3
    return norm_force


def refix(val, p_in="", p_out=""):
    """

    Convert between different SI unit prefixes. Available options are:

    :code:`'T'`   Terra

    :code:`'G'` Giga

    :code:`'M'` Mega

    :code:`'k'` Kilo

    :code:`'m'` Milli

    :code:`'mu'` Micro

    :code:`'n'` Nano

    :code:`'p'` Pico

    Parameters
    ----------
    val: scalar
        The value for which to convert the unit prefix.
    p_in: string, any of the above, optional
        The current prefix of :code:`val`. If :code:`p_in` is undefined,
        :code:`val` has no SI unit prefix.
    p_out: string, any of the above, optional
        The prefix of :code:`val_refix` after the conversion.  If :code:`p_in`
        is undefined, :code:`val_refix` has no SI unit prefix.

    Returns
    -------
    val_refix: scalar
        The value in units of prefix :code:`p_out`.

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
    val_refix = val * prefix[p_in] / prefix[p_out]
    return val_refix


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


if __name__ == "__main__":
    pass
