"""

Test cases for tribology methods

"""

import unittest

import numpy as np

from ..tribology import profball, profrevolve
from .. import tribology_hertz as th
from .. import tribology_lubrication as tl
from .. import tribology_boundary_element as tb


class TestTribology(unittest.TestCase):
    """
    test case methods for general tribology methods
    """

    def test_dyn2kin(self):
        """
        base test for dyn2kin method
        """
        self.assertEqual(tl.dyn2kin(1, 0.5), 0.5)

    def test_kin2dyn(self):
        """
        base test for kin2dyn method
        """
        self.assertEqual(tl.kin2dyn(1, 0.5), 2)

    def test_effective_modulus(self):
        """
        base test for meff method using steel data
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        self.assertEqual(round(e_eff), 230769)

    def test_effective_radii_basic(self):
        """
        base test for eeff method
        """
        r_eff, r_eff_x, r_eff_y = th.reff(6, 0, float('inf'), 3)
        self.assertEqual([r_eff, r_eff_x, r_eff_y], [2, 6, 3])


class TestBoundaryElement(unittest.TestCase):
    """
    test case methods for tribology methods relate to boundary element codes
    """
    def test_be(self):
        """
        base test for combination of boundary element methods
        """
        # inputs for steel ball geometry in contact with steel flat
        r_ball = 6.35
        f_outer = 10
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, _ = th.reff(r_ball, r_ball, 0, 0)

        # create 3d profile for ball
        width, _, _ = th.ahertz(r_eff, r_eff_x, r_eff_x, e_eff, f_outer)
        ax_x, delta_x = np.linspace(-width * 1.1, width * 1.1, 51, retstep=True)
        prof_ball_3d, _ = profrevolve(profball(ax_x, r_ball), ax_x, 2 * r_ball)

        # calculate influence matrix and reduce it
        inf_mat_red = tb.beinflumatred(
            tb.beinflumat(ax_x, ax_x, e_eff))

        # solve for pressure and inner force in ball <-> flat contact
        press, _, f_inner, _ = tb.besolve(prof_ball_3d,
                                          np.zeros((len(ax_x), len(ax_x))),
                                          f_outer, inf_mat_red, delta_x,
                                          delta_x)

        # verify that boundary element solution is equal to analytical solution
        p_max_hertz = 1.5 * th.phertz(r_eff, r_eff_x, r_eff_x, e_eff, f_inner)
        print(p_max_hertz, np.amax(press))
        self.assertEqual(round(np.amax(press)), round(p_max_hertz))


class TestDowsonHamrock(unittest.TestCase):
    """
    test case methods for tribology methods relate to Dowson-Hamrock
    """
    pass


class TestLubrication(unittest.TestCase):
    """
    test case methods for tribology methods relate to Lubrication
    """
    pass


class TestHertz(unittest.TestCase):
    """
    test case methods for tribology methods relate to Hertz contact theory
    """

    def test_hertz_mean_pressure_ball(self):
        """
        base test for phertz
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = th.reff(6.35, 6.35, 0, 0)
        p_mean = th.phertz(r_eff, r_eff_x, r_eff_y, e_eff, 10)
        self.assertEqual(round(p_mean), 574)

    def test_hertz_load_carrying(self):
        """
        base test for fhertz
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = th.reff(6.35, 6.35, 0, 0)
        f_crit = th.fhertz(r_eff, r_eff_x, r_eff_y, e_eff, 574)
        self.assertEqual(round(f_crit, 2), 9.99)

    def test_hertz_half_axes(self):
        """
        base test for ahertz
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = th.reff(15, 0, 0, 10)
        ax_a, ax_b, area = th.ahertz(r_eff, r_eff_x, r_eff_y, e_eff,
                                     100)
        rounded = [round(var, 3) for var in (ax_a, ax_b, area)]
        self.assertEqual(rounded, [0.175, 0.227, 0.125])
