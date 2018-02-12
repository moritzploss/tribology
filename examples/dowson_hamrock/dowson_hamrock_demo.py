"""

Short demonstration of how to use dowson-hamrock function for point contact.

"""

import tribology as tr
import numpy as np

import matplotlib.pyplot as plt


def plot_results(films, film_calc, speeds):
    """

    Plot the calculated and measured data in a loglog plot.

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

    """
    # load data
    database = np.load('demo_ehd_data.npz')
    speeds = database['rolling_speed_mm_s']
    loads = database['ball_load_n']
    films = database['average_film_nm']
    temp = np.mean(database['lube_temp_c'])

    # get material and lube properties
    alpha_p = tr.refix(tr.PressVisc.MINERAL_OIL_GENERIC.value, 'M', '')
    kin = (67.0, 18.9)
    lube_temps = tr.LubeVisc.TEMPS.value
    density = tr.LubeDens.MINERAL_OIL_GENERIC.value

    # calculate dynamic viscosity at operating temperature
    kin_temp = tr.walther(lube_temps[0], kin[0], lube_temps[1], kin[1], temp)
    dyn_temp = tr.refix(tr.kin2dyn(kin_temp, density), 'm', 'M')

    # calculate effective modulus and radius
    e_eff = tr.eeff(tr.YoungsMod.STEEL.value,
                    tr.PoissonRatio.STEEL.value,
                    tr.YoungsMod.GLASS.value,
                    tr.PoissonRatio.GLASS.value)
    r_eff, _, _ = tr.reff(tr.RadBall.TQInch.value, 0, 0, 0)

    # calculate and plot film thickness
    film_calc = tr.dowhampoint(speeds, loads, alpha_p, e_eff, r_eff,
                               dyn_temp)
    plot_results(films, tr.refix(film_calc, 'm', 'n'), speeds)


if __name__ == "__main__":
    dowson_hamrock_demo()
