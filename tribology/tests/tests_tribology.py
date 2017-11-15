"""
Test cases for tribology methods
"""

import unittest

from .. import tribology as tr
from .. import tribology_dowson_hamrock as td
from .. import tribology_hertz as th
from .. import tribology_lubrication as tl


class TestUM(unittest.TestCase):
    """
    test case methods for different tribology methods
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
        base test for effective_modulus method using steel data
        """
        e_eff = tr.effective_modulus(210000, 0.3, 210000, 0.3)
        self.assertEqual(round(e_eff), 230769)

    def test_effective_radii_basic(self):
        """
        base test for effective_radii method
        """
        r_eff, r_eff_x, r_eff_y = tr.effective_radii(6, 0, float('inf'), 3)
        self.assertEqual([r_eff, r_eff_x, r_eff_y], [2, 6, 3])

    def test_hertz_mean_pressure_ball(self):
        """
        base test for hertz_mean_pressure
        """
        e_eff = tr.effective_modulus(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = tr.effective_radii(6.35, 6.35, 0, 0)
        p_mean = th.hertz_mean_pressure(r_eff, r_eff_x, r_eff_y, e_eff, 10)
        self.assertEqual(round(p_mean), 574)

    def test_hertz_load_carrying(self):
        """
        base test for hertz_load_carrying
        """
        e_eff = tr.effective_modulus(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = tr.effective_radii(6.35, 6.35, 0, 0)
        f_crit = th.hertz_load_carrying(r_eff, r_eff_x, r_eff_y, e_eff, 574)
        self.assertEqual(round(f_crit, 2), 9.99)

    def test_hertz_half_axes(self):
        """
        base test for hertz_half_axes
        """
        e_eff = tr.effective_modulus(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = tr.effective_radii(15, 0, 0, 10)
        ax_a, ax_b, area = th.hertz_half_axes(r_eff, r_eff_x, r_eff_y, e_eff,
                                              100)
        rounded = [round(var, 3) for var in (ax_a, ax_b, area)]
        self.assertEqual(rounded, [0.175, 0.227, 0.125])
