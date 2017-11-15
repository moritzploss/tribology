"""
Short demonstration of how to use dowson-hamrock equation for point contact
"""

import matplotlib.pyplot as plt
import numpy as np

from tribology import AlphaP, YoungsMod, PoissonRatio, LubeViscosity, \
    LubeDensity, BallRadii, dowson_hamrock_point, kin2dyn, convert_prefix, \
    walther, effective_modulus, effective_radii


def plot_results(films, film_calc, speeds):
    """
    plot calculated and measured data in loglog plot
    :param films: measured film data, vector, nm
    :param film_calc: calculated film data, vector, nm
    :param speeds: entrainment speeds, mm/s
    :return: None
    """
    plt.figure()
    plt.loglog(speeds, films, 'ko', label='measured', color='black')
    plt.loglog(speeds, film_calc, 'ko', label='calculated', color='grey')
    plt.xlabel('entrainment speed in mm s$^{-1}$')
    plt.ylabel('film thickness in nm')
    plt.legend()
    plt.tight_layout()
    plt.ylim([10, 5000])
    plt.show()


def dowson_hamrock_demo():
    """
    Short demo function to show how to use tribology toolbox for dowson-hamrock
    film thickness calculation. EHD input data file in npz-format required.
    :return: None
    """
    # load data
    database = np.load('demo_ehd_data.npz')
    speeds = database['rolling_speed_mm_s']
    loads = database['ball_load_n']
    films = database['average_film_nm']
    temp = np.mean(database['lube_temp_c'])

    # get material and lube properties
    alpha_p = convert_prefix(AlphaP.MINERAL_OIL_GENERIC.value, 'M', '')
    kin = LubeViscosity.SIGMA_ALDRICH_MINERAL_OIL_HEAVY.value
    lube_temps = LubeViscosity.TEMPS.value
    density = LubeDensity.SIGMA_ALDRICH_MINERAL_OIL_HEAVY.value

    # calculate dynamic viscosity at operating temperature
    kin_temp = walther(lube_temps[0], kin[0], lube_temps[1], kin[1], temp)
    dyn_temp = convert_prefix(kin2dyn(kin_temp, density), 'm', 'M')

    # calculate effective modulus and radius
    e_eff = effective_modulus(YoungsMod.STEEL.value, PoissonRatio.STEEL.value,
                              YoungsMod.GLASS.value, PoissonRatio.GLASS.value)
    r_eff, _, _ = effective_radii(BallRadii.ThreeQuarterInch.value, 0, 0, 0)

    # calculate and plot film thickness
    film_calc = dowson_hamrock_point(speeds, loads, alpha_p, e_eff, r_eff,
                                     dyn_temp)
    plot_results(films, convert_prefix(film_calc, 'm', 'n'), speeds)


if __name__ == "__main__":
    dowson_hamrock_demo()
