# -*- coding: utf-8 -*-

"""

Methods related to Dowson-Hamrock film thickness calculations

"""


def __dowson_hamrock_parameters(r_eff, param_g, param_u, param_w):
    """

    part of dowson-hamrock

    :param r_eff: effective radius
    :param param_g: elasticity_parameter
    :param param_u: velocity_parameter
    :param param_w: load parameter

    :return: ehd parameter

    """
    return r_eff * param_g ** 0.53 * param_u ** 0.67 * param_w ** -0.067


def edowham(alpha_p, e_eff):
    """

    calculate elasticity parameter of dowson-hamrock equation

    :param alpha_p: pressure-viscosity coefficient
    :param e_eff: effective young's modulus

    :return: elasticity parameter

    """
    return alpha_p * e_eff


def vdowham(eta, speed, e_eff, r_eff):
    """

    calculate velocity parameter of dowson-hamrock equation

    :param eta: dynamic viscosity
    :param speed: entrainment speed
    :param e_eff: effective young's modulus
    :param r_eff: effective radius

    :return: velocity parameter

    """
    return eta * speed / (e_eff * r_eff)


def dowhamline(speed, force, alpha_p, e_eff, r_eff, eta, l_eff):
    """

    Calculate mean film thickness according to Dowson-Hamrock equation.

    :param speed: entrainment speed, vector or int/float
    :param force: force, vector or int/float
    :param e_eff: effective modulus
    :param alpha_p: pressure-viscosity coefficient in, int/float
    :param eta: dynamic viscosity of lube, int/float
    :param r_eff: effective radius
    :param l_eff: effective length

    :return: mean film thickness, vector or int/float

    """
    param_g = edowham(alpha_p, e_eff)
    param_u = vdowham(eta, speed, e_eff, r_eff)
    param_w = force / (l_eff * r_eff * e_eff)
    return 2.69 * __dowson_hamrock_parameters(r_eff, param_g, param_u, param_w)


def dowhampoint(speed, force, alpha_p, e_eff, r_eff, eta):
    """

    Calculate mean film thickness according to Dowson-Hamrock equation.

    :param speed: entrainment speed, vector or int/float
    :param force: force, vector or int/float
    :param e_eff: effective modulus in MPa
    :param alpha_p: pressure-viscosity coefficient, int/float
    :param eta: dynamic viscosity of lube, int/float
    :param r_eff: effective radius

    :return: mean film thickness, vector or int/float

    """
    param_g = edowham(alpha_p, e_eff)
    param_u = vdowham(eta, speed, e_eff, r_eff)
    param_w = force / (r_eff ** 2 * e_eff)
    return 1.9 * __dowson_hamrock_parameters(r_eff, param_g, param_u, param_w)
