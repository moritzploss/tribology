import math

import numpy as np

from Constants import SubDir
from system_functions import make_directory
import matplotlib.pyplot as plt


def generate_pub_figures(sys, sim):
    """Generate summary figure for Cylindrical Roller Bearing Tribosystems"""
    plot_axis_pol = np.linspace(-math.pi, math.pi - 2 * math.pi / sys.res_pol,
                                int(sys.res_pol))
    pv_max_ring1 = np.concatenate(
        (sys.pv_ring1[round(sys.res_pol / 2) - 1:-1, :],
         sys.pv_ring2[0:round(sys.res_pol / 2), :]), axis=0)
    eakin_data = np.concatenate(
        (sys.e_akin_ring2.transpose()[round(sys.res_pol / 2) - 1:-1, :],
         sys.e_akin_ring2.transpose()[0:round(sys.res_pol / 2), :]), axis=0)
    pressure_data = np.concatenate(
        (sys.max_pressure_ring1[round(sys.res_pol / 2) - 1:-1, :],
         sys.max_pressure_ring1[0:round(sys.res_pol / 2), :]), axis=0)

    plt.figure()
    ax = plt.subplot(111, polar=True, facecolor='White')
    ax.grid(True)
    ax.plot(plot_axis_pol,
            pv_max_ring1[:, round(sys.roller.res_x / 2) - 1] / np.amax(
                pv_max_ring1), 'k-', label='pv')
    ax.plot(plot_axis_pol,
            eakin_data[:, round(sys.roller.res_x / 2) - 1] / np.amax(
                eakin_data), 'k--', label='eakin')
    ax.plot(plot_axis_pol,
            pressure_data[:, round(sys.roller.res_x / 2) - 1] / np.amax(
                pressure_data), 'k-.',
            label='max pressure')
    ax.plot(sys.pol_ax, sys.normal_forces / np.amax(sys.normal_forces), 'k:',
            label='normal force')
    ax.set_theta_zero_locatioaxisbgn("S")
    ax.set_rlabel_position(112.5)
    label_position = ax.get_rlabel_position()
    r_max = 1.1
    r_label = 'normalised value'
    ax.text(math.radians(label_position + 15), r_max / 2, r_label,
            rotation=22.5, ha='center', va='center')
    ax.set_rmax(r_max)
    make_directory(sim.results_folder, SubDir.pub_figs.value)
    file_handle = sim.results_folder + '/' + SubDir.pub_figs.value + '/' + \
                  'normalised-inner-ring' + '.png'
    plt.savefig(file_handle, bbox_inches='tight', dpi=300)
    plt.close()
