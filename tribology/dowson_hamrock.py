# -*- coding: utf-8 -*-

"""

This module contains functions related to Dowson-Hamrock film thickness
calculations.

"""


def __dowson_hamrock_parameters(r_eff, param_g, param_u, param_w):
    """

     Calculate the EHD-parameter for Dowson-Hamrock film thickness calculations.
     The name EHD-parameter is used within the tribology package to refer to
     what is calculated below; this parameter is not officially defined by
     Dowson-Hamrock.

    Parameters
    ----------
    r_eff: scalar
        The effective radius of the contact problem.
    param_g: scalar
        The elasticity parameter of the contact problem.
    param_u: scalar
        The velocity parameter of the contact problem.
    param_w: scalar
        The load parameter of the contact problem.

    Returns
    -------
    param_ehd: scalar
        The EHD parameter of the contact problem.

    """
    param_ehd = r_eff * param_g ** 0.53 * param_u ** 0.67 * param_w ** -0.067
    return param_ehd


def edowham(alpha_p, e_eff):
    """

    Calculate the elasticity parameter of the contact problem according to
    Dowson-Hamrock.

    Parameters
    ----------
    alpha_p: ndarray, scalar
        The pressure-viscosity coefficient of the lubricant.
    e_eff: ndarray, scalar
        The effective modulus of the contact problem.

    Returns
    -------
    param_elasticity: ndarray, scalar
        The elasticity parameter of the contact problem.

    """
    param_elasticity = alpha_p * e_eff
    return param_elasticity


def vdowham(eta, vel_entrain, e_eff, r_eff):
    """

    Calculate the velocity parameter of the contact problem according to
    Dowson-Hamrock.

    Parameters
    ----------
    eta: ndarray, scalar
        The dynamic viscosity of the lubricant.
    vel_entrain: ndarray, scalar
        The entrainment velocity of the contact problem.
    e_eff: ndarray, scalar
        The effective modulus of the contact problem.
    r_eff: ndarray, scalar
        The effective radius of the contact problem.

    Returns
    -------
    param_velocity: ndarray, scalar
        The velocity parameter of the contact problem.

    """
    param_velocity = eta * vel_entrain / (e_eff * r_eff)
    return param_velocity


def dowhamline(vel_entrain, force, alpha_p, e_eff, r_eff, eta, l_eff):
    """

    Calculate the mean film thickness in a line contact according to
    Dowson-Hamrock.

    Parameters
    ----------
    vel_entrain: ndarray, scalar
        The entrainment velocity of the contact problem.
    force: ndarray, scalar
        The normal force in the contact.
    alpha_p: ndarray, scalar
        The pressure-viscosity coefficient of the lubricant.
    e_eff: ndarray, scalar
        The effective modulus of the contact problem.
    r_eff: ndarray, scalar
        The effective radius of the contact problem.
    eta: ndarray, scalar
        The dynamic viscosity of the lubricant.
    l_eff: ndarray, scalar
        The effective length of the contact.

    Returns
    -------
    h_0: ndarray, scalar
        The central lubricating film thickness in the contact.

    """
    param_g = edowham(alpha_p, e_eff)
    param_u = vdowham(eta, vel_entrain, e_eff, r_eff)
    param_w = force / (l_eff * r_eff * e_eff)
    h_0 = 2.69 * __dowson_hamrock_parameters(r_eff, param_g, param_u, param_w)
    return h_0


def dowhampoint(vel_entrain, force, alpha_p, e_eff, r_eff, eta):
    """

    Calculate the mean film thickness in a point contact according to
    Dowson-Hamrock.

    Parameters
    ----------
    vel_entrain: ndarray, scalar
        The entrainment velocity of the contact problem.
    force: ndarray, scalar
        The normal force in the contact.
    alpha_p: ndarray, scalar
        The pressure-viscosity coefficient of the lubricant.
    e_eff: ndarray, scalar
        The effective modulus of the contact problem.
    r_eff: ndarray, scalar
        The effective radius of the contact problem.
    eta: ndarray, scalar
        The dynamic viscosity of the lubricant.

    Returns
    -------
    h_0: ndarray, scalar
        The central lubricating film thickness in the contact.

    """
    param_g = edowham(alpha_p, e_eff)
    param_u = vdowham(eta, vel_entrain, e_eff, r_eff)
    param_w = force / (r_eff ** 2 * e_eff)
    h_0 = 1.9 * __dowson_hamrock_parameters(r_eff, param_g, param_u, param_w)
    return h_0
