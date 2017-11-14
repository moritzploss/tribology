import unittest
from tribology import dyn2kin, kin2dyn, hertz_mean_pressure, effective_radii, \
    effective_modulus


class TestUM(unittest.TestCase):

    def test_dyn2kin(self):
        self.assertEqual(dyn2kin(1, 0.5), 0.5)

    def test_kin2dyn(self):
        self.assertEqual(kin2dyn(1, 0.5), 2)

    def test_effective_modulus(self):
        e_eff = effective_modulus(210000, 0.3, 210000, 0.3)
        self.assertEqual(round(e_eff), 230769)

    def test_effective_radii_basic(self):
        r_eff, r_eff_x, r_eff_y = effective_radii(6, 0, float('inf'), 3)
        self.assertEqual([r_eff, r_eff_x, r_eff_y], [2, 6, 3])

    def test_hertz_mean_pressure_ball(self):
        e_eff = effective_modulus(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = effective_radii(6.35, 6.35, 0, 0)
        p_mean = hertz_mean_pressure(r_eff, r_eff_x, r_eff_y, e_eff, 10)
        self.assertEqual(round(p_mean), 574)

