"""

Test cases for tribology functions

"""
import unittest
import os
import math
import sys

import numpy as np

sys.path.insert(0, "tribology/p3can")

from ..tribology import profball, profrevolve, profrolleriso
from .. import hertz as th
from .. import lubrication as tl
from .. import boundary_element as tb
from .. import data_import as td
from .. import roller_bearings as trb
from .. import process_slim_mapper as psm
from .. import rough_surfaces as rs
from ..p3can.p3can import p3can


class TestTribology(unittest.TestCase):
    """
    test case methods for general tribology functions
    """

    def test_dyn2kin(self):
        """
        base test for dyn2kin method
        """
        self.assertEqual(tl.dyn2kin(1, 0.5), 2)

    def test_kin2dyn(self):
        """
        base test for kin2dyn method
        """
        self.assertEqual(tl.kin2dyn(1, 0.5), 0.5)

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

    def test_merged_npz(self):
        """
        check if npz database files can be merged
        """
        file_1, _, _ = td.import_del(self.demo_1_txt)
        file_2, _, _ = td.import_del(self.demo_2_txt)
        files = [file_1, file_2]
        merged = td.merge_npz(files)
        self.assertEqual(merged['fx_n'][35], 0.20669)
        for file in files:
            os.remove(file)

    def test_merged_accum_npz(self):
        """
        check if npz accumulated data can be generated during import
        """
        file_1, _, _ = td.import_del(self.demo_1_txt)
        file_2, _, _ = td.import_del(self.demo_2_txt)
        files = [file_1, file_2]
        merged = td.merge_npz(files, accum=['fx_n'])
        self.assertEqual(merged['fx_n'][35], 0.41338)
        for file in files:
            os.remove(file)

    def test_merged_safe_keys_npz(self):
        """
        check that safe merge throws exception if database keys don't
        match
        """
        file_1, _, _ = td.import_del(self.demo_1_txt)
        file_2, _, _ = td.import_del(self.demo_1_pcs_txt)
        files = [file_1, file_2]
        self.assertRaises(KeyError, td.merge_npz, files, accum=['foo'])
        for file in files:
            os.remove(file)

    def test_merged_safe_accum_npz(self):
        """
        check that safe merge throws exception if key(s) in accum are not in
        database
        """
        file_1, _, _ = td.import_del(self.demo_1_txt)
        files = [file_1, file_1]
        self.assertRaises(KeyError, td.merge_npz, files, accum=['foo'])
        os.remove(file_1)

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

    def test_import_multineline_header(self):
        """
        check if multiline column header import works as expected
        """
        f_out, status, _ = td.import_del(
            'tribology/tests/data_import/multicolumn/demo_multiline_header.txt',
            colheadlines=2)
        self.assertEqual(status, True)
        dat = np.load(f_out)
        self.assertEqual(dat['force_n'][2], 8)
        os.remove(f_out)

    def test_import_trios(self):
        """
        check if trios output file can be imported successfully
        """
        f_out, status, _ = td.import_del(
            'tribology/tests/data_import/multicolumn/trios_demo.txt',
            colheadlines=2)
        self.assertEqual(status, True)
        dat = np.load(f_out)
        self.assertEqual(dat['viscosity_pa_s'][50], 0.0363864)
        os.remove(f_out)

    def test_split_del(self):
        """
        check if delimited text file can be split into several files
        """
        outfiles = td.split_del('tribology/tests/data_import/demo_1.txt',
                                hspan=1)
        self.assertEqual(3, len(outfiles))
        for file in outfiles:
            os.remove(file)


class TestHertz(unittest.TestCase):
    """
    test case methods for tribology functions related to Hertz contact theory
    """

    def test_hertz_mean_pressure_circle(self):
        """
        base test for phertz circle
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = th.reff(6.35, 6.35, 0, 0)
        p_mean = th.phertz(r_eff, r_eff_x, r_eff_y, e_eff, 10)
        self.assertEqual(round(p_mean), 574)

    def test_hertz_max_pressure_circle(self):
        """
        base test for phertz circle
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = th.reff(6.35, 6.35, 0, 0)
        p_mean = th.phertz(r_eff, r_eff_x, r_eff_y, e_eff, 10, ret='max')
        self.assertEqual(round(p_mean), 861)

    def test_hertz_pressure_valueerror(self):
        """
        base test for phertz
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = th.reff(6.35, 6.35, 0, 0)
        with self.assertRaises(ValueError):
            th.phertz(r_eff, r_eff_x, r_eff_y, e_eff, 10, ret='foo')

    def test_hertz_mean_pressure_ellipse(self):
        """
        base test for phertz ellipse
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = th.reff(15, 5, 0, 0)
        p_mean = th.phertz(r_eff, r_eff_x, r_eff_y, e_eff, 500)
        self.assertEqual(round(p_mean / 1000, 2), 1.78)

    def test_hertz_mean_pressure_ellipse_2(self):
        """
        base test for phertz ellipse
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = th.reff(0, 0, 5, 15)
        p_mean = th.phertz(r_eff, r_eff_x, r_eff_y, e_eff, 500)
        self.assertEqual(round(p_mean / 1000, 2), 1.78)

    def test_hertz_mean_pressure_ellipse_pos(self):
        """
        base test for phertz ellipse two bodies negative
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = th.reff(15, 5, 20, 15)
        p_mean = th.phertz(r_eff, r_eff_x, r_eff_y, e_eff, 500)
        self.assertEqual(round(p_mean / 1000, 2), 2.33)

    def test_hertz_mean_pressure_ellipse_neg(self):
        """
        base test for phertz ellipse two bodies negative
        """
        e_eff = th.eeff(210000, 0.3, 210000, 0.3)
        r_eff, r_eff_x, r_eff_y = th.reff(15, 5, -20, -15)
        p_mean = th.phertz(r_eff, r_eff_x, r_eff_y, e_eff, 500)
        self.assertEqual(round(p_mean / 1000, 2), 1.06)

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

    def test_approx_hertz_rad_roller(self):
        """
        base test for approx_hertz_rad, using roller with ISO geometry
        """
        axis = np.linspace(-4.99, 4.99, 201)
        prof = profrolleriso(axis, 9, 10)
        rad, _, _ = th.approx_hertz_rad(axis, prof)
        self.assertEqual(round(rad), 2114)

    def test_approx_hertz_rad_half_roller(self):
        """
        assert that same result is obtained for half ISO roller as for complete
        roller
        """
        axis = np.linspace(0, 4.99, 101)
        prof = profrolleriso(axis, 9, 10)
        rad, _, _ = th.approx_hertz_rad(axis, prof)
        self.assertEqual(round(rad), 2114)

    def test_approx_hertz_rad_pos_neg_prof(self):
        """
        assert that same result is obtained for negative profile as for positive
        roller
        """
        axis = np.linspace(0, 4.99, 101)
        prof = profrolleriso(axis, 9, 10)
        rad, _, _ = th.approx_hertz_rad(axis, prof)
        rad_neg, _, _ = th.approx_hertz_rad(axis, -1 * prof)
        self.assertEqual(rad, rad_neg)

    def test_approx_hertz_rad_pos_neg_ax(self):
        """
        assert that same result is obtained for negative profile as for positive
        roller
        """
        axis = np.linspace(0, 4.99, 101)
        prof = profrolleriso(axis, 9, 10)
        rad, _, _ = th.approx_hertz_rad(axis, prof)
        rad_neg, _, _ = th.approx_hertz_rad(-1 * axis, prof)
        self.assertEqual(rad, rad_neg)

    def test_approx_hertz_straight_line(self):
        """
        assert that approx_hertz_rad function returns inf for straight line
        """
        axis = np.linspace(0, 1, 101)
        prof = np.zeros(101)
        rad, _, _ = th.approx_hertz_rad(axis, prof)
        self.assertEqual(rad, float('Inf'))


class TestRoughSurfaces(unittest.TestCase):
    """
    test case methods for tribology functions related to Hertz contact theory
    """

    def test_randsurf(self):
        """
        base test for randsurf method
        """
        ax_x, delta_x = np.linspace(-20, 20, 51, retstep=True)
        ax_y, delta_y = np.linspace(-20, 20, 101, retstep=True)
        heights = rs.randsurf(len(ax_x), len(ax_y), delta_x, delta_y, 0.1, 0.5,
                              0.5)
        self.assertEqual(np.round(np.std(heights), 5), 0.1)
        self.assertLess(abs(np.mean(heights)), 10**-16)


class TestSlimMapper(unittest.TestCase):
    """
    test case methods for tribology functions related to SLIM mapper processing
    """

    def test_slim2thick(self):
        """
        check that processing of image gives specified thickness and rgb values
        """
        spacer, _, _ = td.import_pcs(
            'tribology/tests/process_slim_mapper/demo-3D_SpacerCalibration.txt')

        thick, rgb, _, _, _, _ = psm.slim2thick(
            'tribology/tests/process_slim_mapper/demo-slim-mapper.bmp',
            spacer,
            skip=5,
            crop=0.2,
            aperture={
                'top': 0.667,
                'bottom': 0.667,
            }
        )

        self.assertEqual(np.round(np.nanmean(thick)), 54.0)
        self.assertEqual(([np.round(c, 3) for c in rgb[31][40]]),
                         [0.525, 0.439, 0.357])
        os.remove(spacer)


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
        """
        Test if P3CAN runs successfully with given template
        """
        _ = p3can('tribology/tests/p3can/' +
                  'Template07_BallOn3Plates.py')
