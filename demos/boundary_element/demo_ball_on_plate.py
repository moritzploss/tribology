"""

Short demonstration of how to use methods related to boundary element theory.

"""

import tribology as tr
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_results(ax_x, ax_y, press):
    """

    Plot the 3D pressure distribution.

    """

    # generate plot grid
    x_grid, y_grid = np.meshgrid(ax_x, ax_y, indexing='ij')

    # generate plot
    fig = plt.figure()
    ax_3d = fig.add_subplot(111, projection='3d')
    ax_3d.plot_surface(x_grid, y_grid, press)
    ax_3d.set_xlabel('contact width in mm')
    ax_3d.set_ylabel('contact length in mm')
    ax_3d.set_zlabel('contact pressure in MPa')

    plt.show()


def demo_ball_on_plate():
    """

    Short demo function to show how to use the tribology package for boundary
    element codes. The code below calculates the contact pressure between a
    1/2 inch steel ball in contact with a flat steel plate. The normal force
    is 10 N; the grid size is 51 x 41 in x and y direction respectively.

    """
    # inputs for steel ball geometry in contact with steel flat
    r_ball = tr.RadBall.HInch.value
    f_outer = 10
    e_eff = tr.eeff(tr.YoungsMod.STEEL.value, tr.PoissonRatio.STEEL.value,
                    tr.YoungsMod.STEEL.value, tr.PoissonRatio.STEEL.value)
    r_eff, r_eff_x, r_eff_y = tr.reff(r_ball, r_ball, 0, 0)

    # create 3d profile for ball and flat plate. calculate Hertz contact width
    # first and use result to determine grid size for boundary element solution
    width, _, _ = tr.ahertz(r_eff, r_eff_x, r_eff_y, e_eff, f_outer)
    ax_x, delta_x = np.linspace(-width * 1.1, width * 1.1, 51, retstep=True)
    ax_y, delta_y = np.linspace(-width * 1.1, width * 1.1, 41, retstep=True)
    ball_3d, _ = tr.profrevolve(tr.profball(ax_x, r_ball), ax_y, 2 * r_ball)
    plate_3d = np.zeros((len(ax_x), len(ax_y)))

    # calculate influence matrix and reduce it
    inf_mat = tr.beinflumat(ax_x, ax_y, e_eff)
    inf_mat_red = tr.beinflumatred(inf_mat)

    # solve for pressure and inner force in ball <-> flat contact
    press, _, f_inner, _ = tr.besolve(ball_3d, plate_3d, f_outer,
                                      inf_mat_red, delta_x, delta_y)

    # calculate analytical solution using Hertz equations
    p_max_hertz = 1.5 * tr.phertz(r_eff, r_eff_x, r_eff_x, e_eff, f_inner)

    # show results
    plot_results(ax_x, ax_y, press)
    print('max pressure hertz: {} MPa\n'
          'max pressure boundary element: {} MPa'
          .format(int(p_max_hertz), int(np.amax(press))))


if __name__ == "__main__":
    demo_ball_on_plate()
