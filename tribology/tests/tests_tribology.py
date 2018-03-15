"""

Test cases for tribology functions

"""
import glob
import unittest
import os
import math

import numpy as np

import sys
sys.path.insert(0, "tribology/p3can")

from ..tribology import profball, profrevolve
from .. import hertz as th
from .. import lubrication as tl
from .. import boundary_element as tb
from .. import data_import as td
from .. import roller_bearings as trb
from ..p3can.p3can import p3can


class TestTribology(unittest.TestCase):
    """
    test case methods for general tribology functions
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
    test case methods for tribology functions related to boundary element codes
    """
    def test_be(self):
        """
        base test for combination of boundary element functions
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
    test case methods for tribology functions related to Dowson-Hamrock
    """
    pass


class TestLubrication(unittest.TestCase):
    """
    test case methods for tribology functions related to Lubrication
    """
    pass


class TestRollerBearings(unittest.TestCase):
    """
    test case methods for functions related to roller_bearings
    """

    def test_fcylrolbear(self):
        """
        calculate load distribution in radial roller bearing and compare result
        to reference result
        """
        f_rols, _, _ = trb.fcylrolbear(np.linspace(0, 2 * math.pi, 13),
                                       np.zeros(31),
                                       np.linspace(-5, 5, 31),
                                       3000,
                                       rad_clear=0.0184,
                                       max_dif=0.0005)
        self.assertEqual(round(np.amax(f_rols)), round(967.7))


class TestDataImport(unittest.TestCase):
    """
    test case methods for functions related to data_import
    """
    demo_1 = 'tribology/tests/data_import/demo_1'
    demo_2 = 'tribology/tests/data_import/demo_2'
    demo_1_pcs = 'tribology/tests/data_import/test_data_pcs/demo_1'
    demo_2_pcs = 'tribology/tests/data_import/test_data_pcs/demo_2'

    demo_1_txt = '{}.txt'.format(demo_1)
    demo_2_txt = '{}.txt'.format(demo_2)
    demo_1_pcs_txt = '{}.txt'.format(demo_1_pcs)
    demo_2_pcs_txt = '{}.txt'.format(demo_2_pcs)

    demo_1_npz = '{}.npz'.format(demo_1)
    demo_2_npz = '{}.npz'.format(demo_2)
    demo_1_pcs_npz = '{}.npz'.format(demo_1_pcs)
    demo_2_pcs_npz = '{}.npz'.format(demo_2_pcs)

    demo_dir = 'tribology/tests/data_import'
    demo_dir_pcs = 'tribology/tests/data_import/test_data_pcs'
    demo_dir_rec = 'tribology/tests/data_import/test_data_recursive'

    def test_import_txt_to_npz(self):
        """
        check if simple txt file can be imported to npz format
        """
        f_out, status, _ = td.import_del(self.demo_1_txt)
        self.assertEqual(status, True)
        os.remove(f_out)

    def test_import_txt_database(self):
        """
        check if database is created correctly
        """
        f_out, _, _ = td.import_del(self.demo_1_txt)
        database = np.load(f_out)
        self.assertEqual(database['fx_n'][0], 0.211219)
        os.remove(f_out)

    def test_import_txt_to_mat(self):
        """
        check if simple txt file can be imported to mat format
        """
        f_out, status, _ = td.import_del(self.demo_1_txt, out_ext='mat')
        self.assertEqual(status, True)
        os.remove(f_out)

    def test_import_dir_txt_to_npz(self):
        """
        check if directory containing delimited txt files can be imported to npz
        format
        """
        f_in, f_out, status = td.import_dir(self.demo_dir)
        self.assertEqual(f_in, [self.demo_1_txt, self.demo_2_txt])
        self.assertEqual(f_out[0].endswith(self.demo_1_npz), True)
        self.assertEqual(f_out[1].endswith(self.demo_2_npz), True)
        self.assertEqual(status, [True, True])
        for file in f_out:
            os.remove(file)

    def test_import_pcs_mtm_to_npz(self):
        """
        check if MTM (PCS Instruments) output file can be imported to npz
        """
        f_out, status, _ = td.import_pcs(self.demo_2_pcs_txt, out_ext='npz')
        self.assertEqual(status, True)
        database = np.load(f_out)
        self.assertEqual(database['step_time_s'][0], 23)
        os.remove(f_out)

    def test_import_pcs_ehd_to_npz(self):
        """
        check if EHD2 (PCS Instruments) output file can be imported to npz
        """
        f_out, status, _ = td.import_pcs(self.demo_1_pcs_txt, out_ext='npz')
        self.assertEqual(status, True)
        database = np.load(f_out)
        self.assertEqual(database['film'][0], 130)
        os.remove(f_out)

    def test_import_dir_pcs_to_npz(self):
        """
        check if directory containing pcs files can be imported to npz format
        """
        _, f_out, status = td.import_dir(self.demo_dir_pcs, pcs=True)
        self.assertEqual(status, [True for _ in f_out])
        database = np.load(self.demo_1_pcs_npz)
        self.assertEqual('film_surf' in database, True)
        for file in f_out:
            os.remove(file)

    def test_import_dir_recursive(self):
        """
        check if directory tree can be imported recursively and output files are
        saved together with input files
        """
        _, f_out, status = td.import_dir(self.demo_dir_rec, recursive=True)
        self.assertEqual(status, [True for _ in f_out])
        self.assertEqual(len(f_out), 4)
        for file in f_out:
            os.remove(file)

    def test_import_dir_recursive_out(self):
        """
        check if directory tree can be imported recursively and output files
        can be written to specified output directory
        """
        _, f_out, status = td.import_dir(self.demo_dir_rec, recursive=True,
                                         out_dir=self.demo_dir)
        self.assertEqual(status, [True for _ in f_out])
        self.assertEqual(len(f_out), 4)
        for file in f_out:
            os.remove(file)


class TestHertz(unittest.TestCase):
    """
    test case methods for tribology functions related to Hertz contact theory
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

    def test_approx_hertz_rad(self):
        """
        base test for approx_hertz_rad
        """
        axis = np.asarray([-1, 0, 1])
        prof = np.asarray([0, 1, 0])
        self.assertEqual(round(th.approx_hertz_rad(axis, prof)), 1)
        prof = np.asarray([0, 0, 0])
        self.assertEqual(th.approx_hertz_rad(axis, prof), float('inf'))


class TestP3can(unittest.TestCase):
    """
    integration tests for p3can. load user input templates and check if
    they run without error.
    """
#
#     def test_template01(self):
#         out_dir = p3can('tribology/tests/p3can/' +
#                         'Template01_SingleRowCylindricalRollerBearing.py')
#
#     def test_template03(self):
#         out_dir = p3can('tribology/tests/p3can/' +
#                         'Template03_CylindricalRollerThustBearing.py')
#
#     def test_template05(self):
#         out_dir = p3can('tribology/tests/p3can/' +
#                         'Template05_PinOnDisk.py')
#
#     def test_template06(self):
#         out_dir = p3can('tribology/tests/p3can/' +
#                         'Template06_4Ball.py')
#
    def test_template07(self):
        out_dir = p3can('tribology/tests/p3can/' +
                        'Template07_BallOn3Plates.py')
#
#     def test_template08(self):
#         out_dir = p3can('tribology/tests/p3can/' +
#                         'Template08_RingOnRing.py')