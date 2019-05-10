# -*- coding: utf-8 -*-

"""

This module contains functions related to lubricants and lubrication.

"""
import math
from math import log10 as lg
from math import pi, e

import numpy as np


def kin2dyn(kin, density):
    """

    Convert from kinematic to dynamic viscosity.

    Parameters
    ----------
    kin: ndarray, scalar
        The kinematic viscosity of the lubricant.
    density: ndarray, scalar
        The density of the lubricant.

    Returns
    -------
    dyn: ndarray, scalar
        The dynamic viscosity of the lubricant.

    """
    dyn = kin * density
    return dyn


def dyn2kin(dyn, density):
    """

    Convert from dynamic to kinematic viscosity.

    Parameters
    ----------
    dyn: ndarray, scalar
        The dynamic viscosity of the lubricant.
    density: ndarray, scalar
        The density of the lubricant.

    Returns
    -------
    kin: ndarray, scalar
        The kinematic viscosity of the lubricant.

    """
    kin = dyn / density
    return kin


def __zedwalther(kin):
    """

    Calculate the z-parameter for the Walther equation (ASTM D341).

    Parameters
    ----------
    kin: scalar
        The kinematic viscosity of the lubricant.

    Returns
    -------
    zed: scalar
        The z-parameter.

    """
    zed = kin + 0.7 + 10 ** (-1.47 - 1.84 * kin - 0.51 * kin ** 2)
    return zed


def __nuwalther(zed):
    """

    Calculate the kinematic viscosity for the Walther equation (ASTM D341).

    Parameters
    ----------
    zed: scalar
        The z-parameter of the Walther equation.

    Returns
    -------
    kin: scalar
        The kinematic viscosity.

    """
    kin = (zed - 0.7) - 10 ** (-0.7487 - 3.295 * (zed - 0.7) +
                               0.6119 * (zed - 0.7) ** 2 - 0.3193 *
                               (zed - 0.7) ** 3)
    return kin


def walther(temp_1, nu_1, temp_2, nu_2, temp_3):
    """

    Calculate the kinematic viscosity at temperature `temp_3` based on the
    kinematic viscosities at temperatures `temp_1` and `temp_2`. The
    implementation follows standard ASTM D341.

    Parameters
    ----------
    temp_1: scalar
        The temperature in :math:`^{\\circ}\\text{C}` that corresponds to the
        kinematic viscosity :code:`nu_1`.
    nu_1: scalar
        The kinematic viscosity in cSt at temperature :code:`temp_1`.
    temp_2: scalar
        The temperature in :math:`^{\\circ}\\text{C}` that corresponds to the
        kinematic viscosity :code:`nu_2`.
    nu_2: scalar
        The kinematic viscosity in cSt at temperature :code:`temp_2`.
    temp_3: ndarray
        The temperature in :math:`^{\\circ}\\text{C}` for which to calculate
        the kinematic viscosity.

    Returns
    -------
    nu_3: ndarray
        The kinematic viscosity in cSt at temperature :code:`temp_3`.

    """
    abs_zero = -273.15
    viscs = [nu_1, nu_2]
    thetas = [temp_1 - abs_zero, temp_2 - abs_zero, temp_3 - abs_zero]

    zed = [__zedwalther(nu) for nu in viscs]
    const_a = (lg(lg(zed[0])) - lg(lg(zed[1])) *
               lg(thetas[0]) / lg(thetas[1])) / \
              (1 - lg(thetas[0]) / lg(thetas[1]))
    const_b = (const_a - lg(lg(zed[1]))) / lg(thetas[1])
    nu_3 = __nuwalther(10 ** 10 ** (const_a - const_b * np.log10(thetas[2])))
    return nu_3


def trheomflat(eta, gap_height, omega, r_a):
    """

    Calculate the torque that occurs during a viscosity measurement using a
    plate-on-plate rheometer.

    Parameters
    ----------
    eta: scalar
        The dynamic viscosity of the fluid.
    gap_height: scalar
        The gap height between the stationary and the rotating plate.
    omega: scalar
        The angular velocity of the rotating plate (in radians per second).
    r_a: scalar
        The outer diameter of the rotating plate.

    Returns
    -------
    torque: scalar
        The torque required to rotate the rotating plate.

    """
    torque = eta * pi * omega * r_a ** 4 / (2 * gap_height)
    return torque


def viscrheomflat(torque, gap_height, omega, r_a):
    """

    Calculate the dynamic viscosity of a fluid based on the measurement torque
    in a plate-on-plate rheometer.

    Parameters
    ----------
    torque: scalar
        The torque required to rotate the rotating plate.
    gap_height: scalar
        The gap height between the stationary and the rotating plate.
    omega: scalar
        The angular velocity of the rotating plate (in radians per second).
    r_a: scalar
        The outer diameter of the rotating plate.

    Returns
    -------
    eta: scalar
        The dynamic viscosity of the fluid.

    """
    eta = 2 * gap_height * torque / (pi * omega * r_a ** 4)
    return eta


def trheomcone(eta, alpha, omega, r_a):
    """

    Calculate the torque that occurs during a viscosity measurement using a
    cone-on-plate rheometer.

    Parameters
    ----------
    eta: scalar
        The dynamic viscosity of the fluid.
    alpha: scalar
        The cone angle in radians.
    omega: scalar
        The angular velocity of the rotating plate (in radians per second).
    r_a: scalar
        The outer diameter of the rotating plate.

    Returns
    -------
    torque: scalar
        The torque required to rotate the rotating plate.

    """
    torque = 2 * pi * eta * omega * r_a ** 3 / (3 * alpha)
    return torque


def viscrheomcone(torque, alpha, omega, r_a):
    """

    Calculate the dynamic viscosity of a fluid based on the measurement torque
    in a cone-on-plate rheometer.

    Parameters
    ----------
    torque: scalar
        The torque required to rotate the rotating plate.
    alpha: scalar
        The cone angle in radians.
    omega: scalar
        The angular velocity of the rotating plate (in radians per second).
    r_a: scalar
        The outer diameter of the rotating plate.

    Returns
    -------
    eta: scalar
        The dynamic viscosity of the fluid.

    """
    eta = 3 * torque * alpha / (2 * pi * omega * r_a ** 3)
    return eta


def barus(eta_0, alpha_p, pressure):
    """

    Calculate the dynamic viscosity at a given pressure based on the dynamic
    viscosity at atmospheric pressure and the pressure-viscosity coefficient.

    Parameters
    ----------
    eta_0: float, ndarray
        The dynamic viscosity at atmospheric pressure.
    alpha_p: float, ndarray
        The pressure-viscosity coefficient of the fluid.
    pressure: float, ndarray
        The pressure at which to calculate the dynamic viscosity.

    Returns
    -------
    eta_p: scalar
        The dynamic viscosity at pressure `pressure`.

    """
    eta_p = eta_0 * np.power(e, (alpha_p * pressure))
    return eta_p


