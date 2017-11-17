# -*- coding: utf-8 -*-

"""

Methods related to lubricants and lubrication

"""


from math import pi, e
from math import log10 as lg


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


def __zedwalther(kin_visc):
    """

    calculate z factor for Walther Equation (ASTM D341)

    :param kin_visc: kinematic viscosity

    :return: z factor

    """
    return kin_visc + 0.7 + \
           10 ** (-1.47 - 1.84 * kin_visc - 0.51 * kin_visc ** 2)


def __nuwalther(zed):
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

    zed = [__zedwalther(nu) for nu in viscs]
    const_a = (lg(lg(zed[0])) - lg(lg(zed[1])) *
               lg(thetas[0]) / lg(thetas[1])) / \
              (1 - lg(thetas[0]) / lg(thetas[1]))
    const_b = (const_a - lg(lg(zed[1]))) / lg(thetas[1])
    return __nuwalther(10 ** 10 ** (const_a - const_b * lg(thetas[2])))


def trheomflat(eta, gap_height, omega, r_a):
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


def viscrheomflat(torque, gap_height, omega, r_a):
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


def trheomcone(eta, alpha, omega, r_a):
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


def viscrheomcone(torque, alpha, omega, r_a):
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
