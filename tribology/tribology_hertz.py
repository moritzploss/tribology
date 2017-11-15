# -*- coding: utf-8 -*-

"""
Methods related to Hertz contact theory
"""

from math import sqrt, pi, log


def auxparamshertz(r_eff_x, r_eff_y):
    """
    Calculate Hertz parameters required for contact area/pressure calculations
    :param r_eff_x: effective radius in x direction
    :param r_eff_y: effective radius in y direction
    :return: dimensionless Hertz parameters
    """
    if r_eff_x == 0 or r_eff_y == 0:
        param_lambda = 0
    else:
        param_lambda = min(r_eff_x / r_eff_y, r_eff_y / r_eff_x)
    kappa = 1 / (1 + sqrt(log(16 / param_lambda) / (2 * param_lambda)) -
                 sqrt(log(4)) + 0.16 * log(param_lambda))
    a_ast = kappa * (1 + (2 * (1 - kappa ** 2)) / (pi * kappa ** 2) -
                     0.25 * log(kappa)) ** (1 / 3)
    b_ast = a_ast / kappa
    return a_ast, b_ast, kappa, param_lambda


def dhertz(e_eff, r_x_1, r_y_1, r_x_2, r_y_2, force):
    """
    Elastic displacement in normal direction for arbitrary bodies according to
    hertz contact theory
    :param e_eff: effective young's modulus in MPa
    :param r_x_1: radius body 1 in x direction in mm
    :param r_y_1: radius body 1 in y direction in mm
    :param r_x_2: radius body 2 in x direction in mm
    :param r_y_2: radius body 2 in y direction in mm
    :param force: normal force in N
    :return: combined normal displacement in mm
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
    return param_a * param_b / r_c * (f_2 / f_1)


def ahertz(r_eff, r_eff_x, r_eff_y, e_eff, force):
    """
    Calculate Hertzian contact area half-axis
    :param r_eff: effective radius of contact bodies
    :param r_eff_x: effective radius in x direction
    :param r_eff_y: effective radius in y direction
    :param e_eff: effective modulus
    :param force: normal load
    :return: half axes and total area
    """
    a_ast, b_ast, _, _ = auxparamshertz(r_eff_x, r_eff_y)
    half_axis_a = a_ast * (3 * force * r_eff / e_eff) ** (1 / 3)
    half_axis_b = b_ast * (3 * force * r_eff / e_eff) ** (1 / 3)
    return half_axis_a, half_axis_b, pi * half_axis_a * half_axis_b


def fhertz(r_eff, r_eff_x, r_eff_y, e_eff, p_critical):
    """
    Calculate load carrying capacity of Hertzian contact
    :param r_eff: effective radius of contact bodies
    :param r_eff_x: effective radius in x direction
    :param r_eff_y: effective radius in y direction
    :param e_eff: effective modulus
    :param p_critical: critical mean contact pressure
    :return: critical normal force
    """
    a_ast, b_ast, _, _ = auxparamshertz(r_eff_x, r_eff_y)
    return (pi * a_ast * b_ast * p_critical) ** 3 * (3 * r_eff / e_eff) ** 2


def phertz(r_eff, r_eff_x, r_eff_y, e_eff, force):
    """
    Calculate mean contact pressure in Hertzian contact
    :param r_eff: effective radius of contact bodies
    :param r_eff_x: effective radius in x direction
    :param r_eff_y: effective radius in y direction
    :param e_eff: effective modulus
    :param force: normal force
    :return: mean contact pressure
    """
    _, _, area = ahertz(r_eff, r_eff_x, r_eff_y, e_eff, force)
    return force / area
