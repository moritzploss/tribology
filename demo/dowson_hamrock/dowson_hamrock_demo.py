"""
Short demonstration of how to use dowson-hamrock equation for point contact
"""

import matplotlib.pyplot as plt
import numpy as np
import tribology as tr


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
    alpha_p = tr.convert_prefix(tr.AlphaP.MINERAL_OIL_GENERIC.value, 'M', '')
    kin = tr.LubeViscosity.SIGMA_ALDRICH_MINERAL_OIL_HEAVY.value
    lube_temps = tr.LubeViscosity.TEMPS.value
    density = tr.LubeDensity.SIGMA_ALDRICH_MINERAL_OIL_HEAVY.value

    # calculate dynamic viscosity at operating temperature
    kin_temp = tr.walther(lube_temps[0], kin[0], lube_temps[1], kin[1], temp)
    dyn_temp = tr.convert_prefix(tr.kin2dyn(kin_temp, density), 'm', 'M')

    # calculate effective modulus and radius
    e_eff = tr.effective_modulus(tr.YoungsMod.STEEL.value,
                                 tr.PoissonRatio.STEEL.value,
                                 tr.YoungsMod.GLASS.value,
                                 tr.PoissonRatio.GLASS.value)
    r_eff, _, _ = tr.effective_radii(tr.RadBall.TQInch.value, 0, 0, 0)

    # calculate and plot film thickness
    film_calc = tr.dowson_hamrock_point(speeds, loads, alpha_p, e_eff, r_eff,
                                        dyn_temp)
    plot_results(films, tr.convert_prefix(film_calc, 'm', 'n'), speeds)


if __name__ == "__main__":
    dowson_hamrock_demo()
