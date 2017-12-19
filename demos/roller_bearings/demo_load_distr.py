"""

Short demonstration of how to calculate the load distribution in a cylindrical
roller bearing using the tribology package.

"""

from math import pi

import matplotlib.pyplot as plt
import numpy as np

import tribology as tr


def load_distr_cylrolbear():
    """

    Calculate the normal force for each roller in a cylindrical roller bearing.
    The example below uses a bearing with 13 rollers and a roller length of 10.
    The rollers and bearing rings have no profile (see comb_prof).

    """

    # calculation inputs
    ang_pos = np.linspace(0, 2 * pi, 13)
    ax_roll = np.linspace(-5, 5, 31)
    rad_clear = 0.0184
    f_rad = 3000
    comb_prof = np.zeros(len(ax_roll))

    # calculate normal force per roller
    f_rolls, _, _ = tr.fcylrolbear(ang_pos, comb_prof, ax_roll, f_rad,
                                   rad_clear=rad_clear, max_dif=0.0005)

    # plot normal force per roller over roller position
    axis = plt.subplot(111, projection='polar')
    axis.plot(ang_pos, f_rolls)
    plt.show()


if __name__ == "__main__":
    load_distr_cylrolbear()
