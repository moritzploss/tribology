import math

import numpy as np

from numeric_methods import secant_method


def load_distr_cyl_rol_bear(phis, number_rollers, length_roller, res_x,
                            combined_profile, roller_axis, radial_clearance,
                            radial_load):
    """Caclulate the load distribution in a cylindrical roller bearing
    according to standard DIN 26281"""
    cl = 35948 * math.pow(length_roller, (8 / 9))
    cs = cl / res_x
    x_axis = abs(roller_axis)
    cos_phi = list((math.cos(phi) for phi in phis))

    delta_r = radial_clearance * 4 + 0.000000001
    delta_f = 20
    sum_re_ns = np.zeros(number_rollers)
    psi_j = np.zeros(number_rollers)
    delta_j = np.zeros(number_rollers)
    sum_mns = np.zeros(number_rollers)
    delta_jk = np.zeros((number_rollers, res_x))
    delta_re = np.zeros(number_rollers)
    x_value = []
    fx_value = []
    radial_load = radial_load or 0.000000001

    while abs(delta_f) > 0.0005 * radial_load:
        sum_zwk = 0
        delta_re = np.zeros(number_rollers)
        for z in range(number_rollers):
            sum_re_ns[z] = 0
            sum_mns[z] = 0
            delta_j[z] = delta_r * cos_phi[z] - radial_clearance / 2
            for n in range(0, res_x):
                delta_jk[z, n] = max((delta_j[z] - x_axis[n] * math.tan(
                    psi_j[z]) - 2 * combined_profile[n]), 0)
                sum_re_ns[z] += math.pow(delta_jk[z, n], (10 / 9))
            delta_re[z] = cos_phi[z] * sum_re_ns[z]
            sum_zwk += delta_re[z]
        delta_f = abs(radial_load - cs * sum_zwk)
        x_value = np.append(x_value, [delta_r])
        fx_value = np.append(fx_value, [delta_f])
        delta_r = secant_method(x_value, fx_value)

    roller_normal_forces = sum_re_ns * cs
    return roller_normal_forces, delta_re
