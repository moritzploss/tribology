import copy
import math
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly
import scipy.io

from BluPrintRollBear import RadialRollBear
from Constants import SubDir, TexTempl, PltOpts, PrintOpts, UnitTex, PreSol
from Constants import Unit
from cartesian_plot_functions import plt_contact, plt_3d, plt_2d_scatt_line, \
    plt_profile_approx, plt_profiles
from generate_latex_output import get_calc_specific_latex_template
from hertz_equations import hertz_displ
from influ_matrix_management import cache_influ_mat, load_influ_mat
from load_distr_cyl_rol_bear import load_distr_cyl_rol_bear
from polar_plot_functions import polplt_scatt_line, polplt_line_line, \
    polplt_line
from solve_half_space import solve_half_space, pre_solve_half_space
from system_functions import make_directory, print_it, print_progress, \
    to_preci, save_to_matlab


class CylRollBear(RadialRollBear):
    """Global bearing data"""

    def __init__(self, number_rollers, roller, ring1, ring2,
                 init_force, radial_clearance=0, res_pol=1,
                 path_roller_slip=None, bearing_name='bearing'):
        super().__init__(bearing_name, number_rollers, roller, ring1, ring2,
                         init_force, radial_clearance, res_pol,
                         path_roller_slip)

        self.ring1.press_zone_len = None
        self.ring2.press_zone_len = None
        self.ring2.eff_footpr_vel = None
        self.ring1.press_pol = None
        self.ring2.press_pol = None
        self.ring2.overroll_t = None
        self.ring2.no_overroll_t = None
        self.ring1.uniq_press = None
        self.ring2.uniq_press = None
        self.slip_data = None

        self.influ_mat_db_1 = None
        self.influ_mat_db_2 = None
        self.roller_norm_forces = None
        self.roller_norm_displ = None
        self.normal_forces = None
        self.norm_displ = None
        self.uniq_norm_displ = None
        self.uniq_norm_forces = None
        self.uniq_force_idx = None

        self.eff_vel_cage = None
        self.pol_coords_slip = None

        self.calc_polar_coordinates()

    def calc_polar_coordinates(self):
        """Set-up a polar coordinate system for the inner/outer ring"""
        if (self.res_pol % self.num_rollers) is not 0:
            self.res_pol = int(
                math.ceil(self.res_pol / self.num_rollers) * self.num_rollers)
            print_it("polar resolution ('res_pol') changed to {}".format(
                self.res_pol), PrintOpts.lvl1.value)

        self.pol_ax = np.linspace(0, 2 * math.pi - 2 * math.pi / self.res_pol,
                                  self.res_pol)
        self.num_roller_pos = int(self.res_pol / self.num_rollers)
        self.d_phi = 2 * math.pi / self.num_rollers
        self.dd_phi = 2 * math.pi / self.res_pol

        self.phi_mat = np.zeros((self.num_roller_pos, self.num_rollers))
        for i in range(self.num_roller_pos):
            self.phi_mat[i, :] = np.linspace(i * self.dd_phi,
                                             2 * math.pi - self.d_phi + i * self.dd_phi,
                                             self.num_rollers)

    def calc_load_distribution(self, ui=None, res_dir=None):
        """Calculate the load distribution for all roller positions by calling
        the function 'load_distr_cyl_rol_bear' while changing the polar
        position of the rollers with respect to the inner/outer ring"""
        print_it(
            'calculating {} load distributions'.format(self.num_roller_pos))

        # initialise variables
        count = 0
        x_profile = abs(self.roller.x_profile + self.ring1.x_profile)
        self.roller_norm_forces = np.zeros(
            (self.num_roller_pos, self.num_rollers))
        self.roller_norm_displ = np.zeros(
            (self.num_roller_pos, self.num_rollers))

        # get all load distributions for all roller positions
        for roller_pos in range(self.num_roller_pos):
            self.roller_norm_forces[roller_pos, :], self.roller_norm_displ[
                                                    roller_pos, :] = \
                load_distr_cyl_rol_bear(self.phi_mat[roller_pos, :],
                                        self.num_rollers, self.roller.length,
                                        self.roller.res_x, x_profile,
                                        self.roller.x_axis, self.rad_clear,
                                        self.init_force)
            count = print_progress(count, self.num_roller_pos)
        self.normal_forces = np.transpose(self.roller_norm_forces).reshape(
            (self.roller_norm_forces.size, 1))
        self.norm_displ = np.transpose(self.roller_norm_displ).reshape(
            (self.roller_norm_displ.size, 1))

        # get unique normal forces and unique normal displacements
        self.uniq_norm_forces, self.uniq_force_idx = np.unique(
            np.round(self.roller_norm_forces), return_inverse=True)
        uniq_norm_displ = np.zeros(len(self.uniq_norm_forces))
        for i in range(len(self.uniq_norm_forces)):
            uniq_norm_displ[i] = np.mean(
                self.norm_displ[self.uniq_force_idx == i])
        self.uniq_norm_displ = sorted(np.absolute(uniq_norm_displ) / 10)

        # save load distribution data to matlab database
        dat_dict = {'roller_positions': self.num_roller_pos,
                    'rollers': self.num_rollers,
                    'roller_normal_forces': self.roller_norm_forces,
                    'uniq_norm_forces': self.uniq_norm_forces,
                    'normal_forces': self.normal_forces,
                    'polar_coordinates': self.phi_mat}
        save_to_matlab(dat_dict, res_dir, 'load-distribution')

    def calc_contact_pressure(self, ui=None, res_dir=None):
        """For each unique normal roller force, calculate the contact pressure
        between inner ring and roller and outer ring and roller"""

        # initialise variables
        print_it('calculating pressure distributions')
        count = 0
        num_uniq_forces = self.uniq_norm_forces.shape[0]
        self.ring1.uniq_press = np.zeros(
            (self.roller.res_x, self.roller.res_y, num_uniq_forces))
        self.ring2.uniq_press = np.zeros(
            (self.roller.res_x, self.roller.res_y, num_uniq_forces))

        # iterate over unique normal forces, starting with highest normal force
        # at end of array
        for unique_norm_force_idx in range(num_uniq_forces - 1, -1, -1):
            if self.uniq_norm_forces[unique_norm_force_idx] > 0:
                norm_force = self.uniq_norm_forces[unique_norm_force_idx]
                norm_displ = self.uniq_norm_displ[unique_norm_force_idx]
                # start with pressure distribution for highest normal force
                # since this determines grid size and
                # influence matrix for all other calculations
                if unique_norm_force_idx == num_uniq_forces - 1:
                    self.solve_highest_pressure(ui, res_dir, norm_force)
                    print_it('solving for {} pressure distributions'.format(
                        2 * num_uniq_forces), PrintOpts.lvl1.value)
                # solve for pressure distributions for other normal forces,
                # reuse influence matrix
                else:
                    self.ring1.loc_press, self.influ_mat_db_1 = \
                        solve_half_space(self.roller.profile,
                                         self.ring1.profile, self.roller.x_axis,
                                         self.roller.y_axis, self.roller.res_x,
                                         self.roller.res_y, self.roller.delta_x,
                                         self.roller.delta_y, self.roller.e,
                                         self.ring1.e, self.roller.ny,
                                         self.ring1.ny, norm_force, res_dir,
                                         influ_mat_db=self.influ_mat_db_1,
                                         print_prog=False,
                                         init_displ=norm_displ)

                    self.ring2.loc_press, self.influ_mat_db_2 = \
                        solve_half_space(self.roller.profile,
                                         self.ring2.profile, self.roller.x_axis,
                                         self.roller.y_axis, self.roller.res_x,
                                         self.roller.res_y, self.roller.delta_x,
                                         self.roller.delta_y, self.roller.e,
                                         self.ring2.e, self.roller.ny,
                                         self.ring2.ny, norm_force, res_dir,
                                         influ_mat_db=self.influ_mat_db_2,
                                         print_prog=False,
                                         init_displ=norm_displ)
                # store pressure data
                self.ring1.uniq_press[:, :,
                unique_norm_force_idx] = self.ring1.loc_press
                self.ring2.uniq_press[:, :,
                unique_norm_force_idx] = self.ring2.loc_press
            count = print_progress(print_progress(count, 2 * num_uniq_forces),
                                   2 * num_uniq_forces)
        self.reassign_pressure_results()
        self.save_pressure_to_matlab(res_dir)
        cache_influ_mat(ui, [self.influ_mat_db_1, self.influ_mat_db_2], res_dir)

    def get_grid_size(self, ui, res_dir, uniq_norm_force):
        """Determine grid size by running (quick) simulation with simplified
        contact bodies. Calculate pressure only for contact between outer ring
        and roller since the contact area will be larger than for the case of
        inner ring and roller."""
        print_it('determining grid size', PrintOpts.lvl1.value)
        self.roller.simple_clone()
        self.roller.clone.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                                       self.init_force)
        self.ring2.simple_clone()
        self.ring2.clone.make_slave_to(self.roller.clone)

        init_displ = hertz_displ(self.roller.clone.e, self.ring2.clone.e,
                                 self.roller.clone.ny, self.ring2.clone.ny,
                                 self.roller.clone.r_hertz_x,
                                 self.roller.clone.r_hertz_y,
                                 self.ring2.clone.r_hertz_x,
                                 self.ring2.clone.r_hertz_y, uniq_norm_force)
        too_many_elements_in_contact = 1
        contact_width = 0.05
        while too_many_elements_in_contact != 0:
            self.roller.clone.make_profile(self.roller.clone.res_x,
                                           self.roller.clone.res_y,
                                           self.init_force,
                                           contact_width=contact_width)
            self.ring2.clone.make_slave_to(self.roller.clone)

            pressure, init_displ = \
                pre_solve_half_space(self.roller.clone.profile,
                                     self.ring2.clone.profile,
                                     self.roller.clone.x_axis,
                                     self.roller.clone.y_axis,
                                     self.roller.clone.res_x,
                                     self.roller.clone.res_y,
                                     self.roller.clone.delta_x,
                                     self.roller.clone.delta_y,
                                     self.roller.clone.e,
                                     self.ring2.clone.e, self.roller.clone.ny,
                                     self.ring2.clone.ny, uniq_norm_force,
                                     init_displ=init_displ, print_prog=False)

            pressure_elements = sum(
                pressure[math.floor(self.roller.clone.res_y / 2), :] > 0)
            too_many_elements_in_contact = self.roller.clone.res_y - \
                                           pressure_elements - 2
            contact_width += -np.sign(
                too_many_elements_in_contact) * contact_width / 25

        self.roller.make_profile(self.roller.res_x, self.roller.res_y,
                                 self.init_force, contact_width=contact_width)
        self.ring1.make_slave_to(self.roller)
        self.ring2.make_slave_to(self.roller)
        return init_displ

    def solve_highest_pressure(self, ui, res_dir, uniq_norm_force):
        """Solve for pressure distribution for highest normal force"""
        init_displ = self.get_grid_size(ui, res_dir, uniq_norm_force)
        [self.influ_mat_db_1, self.influ_mat_db_2] = load_influ_mat(ui, res_dir,
                                                                    2)

        self.ring1.loc_press, self.influ_mat_db_1 = \
            solve_half_space(self.roller.profile, self.ring1.profile,
                             self.roller.x_axis, self.roller.y_axis,
                             self.roller.res_x, self.roller.res_y,
                             self.roller.delta_x, self.roller.delta_y,
                             self.roller.e, self.ring1.e, self.roller.ny,
                             self.ring1.ny, uniq_norm_force, res_dir,
                             print_prog=False, influ_mat_db=self.influ_mat_db_1,
                             init_displ=init_displ)

        self.ring2.loc_press, self.influ_mat_db_2 = \
            solve_half_space(self.roller.profile, self.ring2.profile,
                             self.roller.x_axis, self.roller.y_axis,
                             self.roller.res_x, self.roller.res_y,
                             self.roller.delta_x, self.roller.delta_y,
                             self.roller.e, self.ring2.e, self.roller.ny,
                             self.ring2.ny, uniq_norm_force, res_dir,
                             print_prog=False, influ_mat_db=self.influ_mat_db_2,
                             init_displ=init_displ)

    def reassign_pressure_results(self):
        """Extract and reassign pressure data to allow for easier storage and
        access"""
        idx_mat = np.reshape(self.uniq_force_idx,
                             (self.num_roller_pos, self.num_rollers))
        self.ring1.press_pol = np.zeros(
            (self.roller.res_x, self.roller.res_y, self.res_pol))
        self.ring2.press_pol = copy.copy(self.ring1.press_pol)
        self.ring1.press = np.zeros((self.roller.res_x, self.roller.res_y,
                                     self.num_roller_pos, self.num_rollers))
        self.ring2.press = np.zeros((self.roller.res_x, self.roller.res_y,
                                     self.num_roller_pos, self.num_rollers))
        self.ring1.max_press = np.zeros((self.res_pol, self.roller.res_x))
        self.ring2.max_press = np.zeros((self.res_pol, self.roller.res_x))
        loop_count = 0
        for j in range(self.num_rollers):
            for i in range(self.num_roller_pos):
                self.ring1.press[:, :, i, j] = self.ring1.uniq_press[:, :,
                                               idx_mat[i, j]]
                self.ring2.press[:, :, i, j] = self.ring2.uniq_press[:, :,
                                               idx_mat[i, j]]
                self.ring1.max_press[loop_count, :] = np.amax(
                    self.ring1.press[:, :, i, j], axis=1)
                self.ring2.max_press[loop_count, :] = np.amax(
                    self.ring2.press[:, :, i, j], axis=1)
                self.ring1.press_pol[:, :, loop_count] = self.ring1.uniq_press[
                                                         :, :, idx_mat[i, j]]
                self.ring2.press_pol[:, :, loop_count] = self.ring2.uniq_press[
                                                         :, :, idx_mat[i, j]]
                loop_count += 1

    def calc_kinematics(self, rot_vel1, rot_vel2, ui=None, res_dir=None):
        """Calculate bearing kinematics based on rotational velocities of
        contact bodies"""
        print_it("calculating kinematics")
        self.rot_vel1 = rot_vel1
        self.rot_vel2 = rot_vel2
        self.ring1.omega = rot_vel1 / 60
        self.ring2.omega = rot_vel2 / 60
        self.eff_rot_vel = self.rot_vel1 - self.rot_vel2
        self.eff_omega = self.eff_rot_vel / 60
        self.eff_omega_cage = self.eff_omega / 2

        self.ring1.eff_raceway_vel = self.ring1.diameter * math.pi * \
                                     self.eff_omega
        self.ring2.eff_raceway_vel = 0
        self.eff_vel_cage = (
                            self.ring1.eff_raceway_vel +
                            self.ring2.eff_raceway_vel) / 2
        self.roller.eff_raceway_vel = \
            math.pi * self.roller.diameter * (
            self.mean_diameter / self.roller.diameter) * self.eff_omega / 2
        self.roller.actual_raceway_vel = np.zeros(
            (self.res_pol, self.roller.res_x))
        self.roller.rel_vel = np.zeros((self.res_pol, self.roller.res_x))

        if self.path_roller_slip is None:
            self.roller.actual_raceway_vel[:, :] = self.roller.eff_raceway_vel
        else:
            self.calc_slip_data()

        self.calc_overroll_times()

    def calc_slip_data(self):
        """Load, fit and interpolate splip data from user input file"""
        self.slip_data = []
        with open(self.path_roller_slip, 'r') as f:
            for line in f:
                self.slip_data.append(
                    [float(data_point) for data_point in line.split()])
        self.pol_coords_slip = np.linspace(0, 2 * math.pi, len(self.slip_data))
        slip_fit = poly.polyfit(self.pol_coords_slip, self.slip_data, 17)
        self.roller.slip = poly.polyval(self.pol_ax, slip_fit)[0, :]
        self.roller.slip[np.where(self.roller.slip > 1)[0]] = 1

        for i in range(len(self.pol_ax)):
            self.roller.actual_raceway_vel[i, :] = \
                self.roller.eff_raceway_vel * self.roller.slip[i]
            self.roller.rel_vel[i, :] = \
                self.roller.eff_raceway_vel -\
                self.roller.actual_raceway_vel[i, :]

    def calc_overroll_times(self):
        """Calculate the duration of each overrolling cycle for roller,
        inner ring and outer ring"""
        polar_points_under_load = -1 * np.ones(self.res_pol)
        for i in range(self.res_pol):
            if sum(sum(self.ring2.press_pol[:, :, i])) != 0:
                polar_points_under_load[i] = 1
        circ_ratio_loaded = sum(polar_points_under_load == 1) / self.res_pol
        circ_ratio_unloaded = sum(polar_points_under_load == -1) / self.res_pol

        self.ring1.eff_footpr_vel = 2 * math.pi * (
            self.ring1.diameter / 2) * self.eff_omega_cage
        self.ring2.eff_footpr_vel = 2 * math.pi * (
            self.ring2.diameter / 2) * self.eff_omega_cage
        self.ring1.press_zone_len = (self.ring1.press_pol > 0).sum(
            1) * self.ring1.delta_y
        self.ring2.press_zone_len = (self.ring2.press_pol > 0).sum(
            1) * self.ring2.delta_y

        # case inner ring: if inner ring is standing still
        if self.rot_vel1 == 0:
            self.ring1.overroll_t_incr = \
                self.roller.delta_y / self.ring1.eff_footpr_vel
            self.ring1.no_overroll_t = np.divide(
                ((2 * math.pi * (self.ring1.diameter / 2)) / self.num_rollers -
                 self.ring1.press_zone_len), self.ring1.eff_footpr_vel)
        # case inner ring: if inner ring is standing still
        else:
            t_ratio_loaded_ring1 = circ_ratio_loaded
            t_ratio_unloaded_ring1 = circ_ratio_unloaded
            self.ring1.overroll_t_incr = \
                self.roller.delta_y / self.ring1.eff_footpr_vel
            self.ring1.no_overroll_t = np.divide(
                ((2 * math.pi * (
                self.ring1.diameter / 2)) / self.num_rollers -
                self.ring1.press_zone_len),
                self.ring1.eff_footpr_vel) * t_ratio_loaded_ring1 + 2 * \
                math.pi * (self.ring1.diameter / 2) * t_ratio_unloaded_ring1 / \
                self.ring1.eff_footpr_vel

        # case outer ring: if outer ring is standing still
        if self.rot_vel2 == 0:
            self.ring2.overroll_t_incr = \
                self.roller.delta_y / self.ring2.eff_footpr_vel
            self.ring2.no_overroll_t = np.divide(
                ((2 * math.pi * (self.ring2.diameter / 2)) / self.num_rollers -
                 self.ring2.press_zone_len), self.ring2.eff_footpr_vel)
        # case outer ring: if outer ring is rotating
        else:
            t_ratio_loaded_ring2 = circ_ratio_loaded
            t_ratio_unloaded_ring2 = circ_ratio_unloaded
            self.ring2.overroll_t_incr = \
                self.roller.delta_y / self.ring2.eff_footpr_vel
            self.ring2.no_overroll_t = np.divide(
                ((2 * math.pi * (self.ring2.diameter / 2)) / self.num_rollers -
                self.ring2.press_zone_len),
                self.ring2.eff_footpr_vel) * t_ratio_loaded_ring2 + \
                2 * math.pi * (
                self.ring2.diameter / 2) * t_ratio_unloaded_ring2 / \
                self.ring2.eff_footpr_vel

    def calc_pv(self, ui=None, res_dir=None):
        """Calculate product of local maximum pressure and local maximum
        relative velocity"""
        print_it("calculating pv_rel")
        self.ring1.pv = np.absolute(
            np.multiply(self.ring1.max_press, self.roller.rel_vel) / 1000)
        self.ring2.pv = np.absolute(
            np.multiply(self.ring2.max_press, self.roller.rel_vel) / 1000)
        dat_dict = dict(x_axis=self.roller.x_axis,
                        pv_rel_ir_3d=self.ring1.pv,
                        pv_rel_ir_2d=np.amax(self.ring1.pv, axis=1),
                        pv_rel_or_3d=self.ring2.pv,
                        pv_rel_or_2d=np.amax(self.ring2.pv, axis=1),
                        axis_pol=self.pol_ax)
        save_to_matlab(dat_dict, res_dir, 'pv_rel_ring2')

    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation in W per m^2"""
        print_it("calculating e_a,kin")

        pv_loc_ring1 = np.multiply(self.ring1.press_pol.sum(1),
                                   self.roller.rel_vel.transpose())
        pv_loc_ring2 = np.multiply(self.ring2.press_pol.sum(1),
                                   self.roller.rel_vel.transpose())

        if self.rot_vel2 == 0:
            self.ring2.e_akin = np.absolute(
                np.divide(np.multiply(pv_loc_ring2, self.ring2.overroll_t_incr),
                          self.ring2.no_overroll_t)) / 1000
        else:
            self.ring2.e_akin = np.absolute(np.divide(
                np.multiply(pv_loc_ring2.sum(1), self.ring2.overroll_t_incr),
                self.ring2.no_overroll_t.sum(1))) / 1000
        if self.rot_vel1 == 0:
            self.ring1.e_akin = np.absolute(
                np.divide(np.multiply(pv_loc_ring1, self.ring2.overroll_t_incr),
                          self.ring1.no_overroll_t)) / 1000
        else:
            self.ring1.e_akin = np.absolute(np.divide(
                np.multiply(pv_loc_ring1.sum(1), self.ring1.overroll_t_incr),
                self.ring1.no_overroll_t.sum(1))) / 1000

    def plot_it(self, ui=None, res_dir=None):
        """Orchestrate output plot generation"""
        print_it("plotting results and finishing up")
        plt_profiles(res_dir)
        self.plt_profiles_contacts(res_dir)
        self.plt_pol_slip(res_dir)
        self.plt_load_distr(res_dir)
        self.plt_pres(res_dir)
        self.plt_pv(res_dir)
        self.plt_eakin(res_dir)
        self.generate_summary_figure(res_dir)

    def save_pressure_to_matlab(self, res_dir):
        dat_dict = dict(contact_pressures_IR=self.ring1.press,
                        contact_pressures_OR=self.ring2.press,
                        x_axis=self.roller.x_axis,
                        y_axis=self.roller.y_axis,
                        roller_number=np.linspace(1, self.num_rollers,
                                                  self.num_rollers),
                        roller_position=np.linspace(1, self.num_roller_pos,
                                                    self.num_roller_pos))
        save_to_matlab(dat_dict, res_dir, 'pressure-field')

    def generate_latex_output(self, calc_spec_tex_file_handle, sim, ui=None,
                              res_dir=None):
        """Generate the calculation-specific part of the latex output file"""
        av_press_ring1 = np.mean(
            self.ring1.press[:, :, 0, 0][self.ring1.press[:, :, 0, 0] > 0])
        av_press_ring2 = np.mean(
            self.ring2.press[:, :, 0, 0][self.ring2.press[:, :, 0, 0] > 0])
        out_data = [
            ('pv_rel inner ring, max.', to_preci(np.amax(self.ring1.pv), 4),
             UnitTex.pvrel.value, 'unverified'),
            ('pv_rel outer ring, max.', to_preci(np.amax(self.ring2.pv), 4),
             UnitTex.pvrel.value, 'unverified'),
            ('e_akin inner ring, max.',
             to_preci(np.amax([self.ring1.e_akin]), 4),
             UnitTex.eakin.value, 'unverified'),
            ('e_akin outer ring, max.',
             to_preci(np.amax([self.ring2.e_akin]), 4),
             UnitTex.eakin.value, 'unverified'),
            ('pressure inner ring, max. av.', to_preci(av_press_ring1, 4),
             UnitTex.pressure.value, 'unverified'),
            ('pressure outer ring, max. av.', to_preci(av_press_ring2, 4),
             UnitTex.pressure.value, 'unverified'),
            (
            'normal force, max.', to_preci(np.amax(self.roller_norm_forces), 4),
            UnitTex.force.value, 'unverified'),
            ('pressure inner ring, max. max.',
             to_preci(np.amax(self.ring1.press[:, :, 0, 0]), 4),
             UnitTex.pressure.value, 'unverified'),
            ('pressure outer ring, max. max.',
             to_preci(np.amax(self.ring2.press[:, :, 0, 0]), 4),
             UnitTex.pressure.value, 'unverified'),
            ('contact area inner ring, max.',
             to_preci(self.roller.get_area(self.ring1.press[:, :, 0, 0]), 4),
             UnitTex.area.value, 'unverified'),
            ('contact area outer ring, max.',
             to_preci(self.roller.get_area(self.ring2.press[:, :, 0, 0]), 4),
             UnitTex.area.value, 'unverified')
            ]

        scale_factor_ir = 0.5 + 0.1 * (self.rot_vel1 != 0)
        scale_factor_or = 0.5 + 0.1 * (self.rot_vel2 != 0)

        table_calc_summary = []
        for key, value, unit, status in sorted(out_data):
            table_calc_summary.append(
                (re.sub('_', '\_', key), value, unit, status))

        latex_variables = {'contact_plot1': '{}{}contact1.png'.format(
            SubDir.tex_figs_rel_to_tex_file.value, '/'),
                           'contact_plot2': '{}{}contact2.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/'),
                           'scale_factor_ir': scale_factor_ir,
                           'scale_factor_or': scale_factor_or,
                           'table_calc_summary': table_calc_summary,
                           'load_plot1': '{}{}load1.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/'),
                           'pressure_plot1': '{}{}pressure1.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/'),
                           'pressure_plot2': '{}{}pressure2.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/'),
                           'slip_plot1': '{}{}slip1.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/'),
                           'energy_plot1': '{}{}energy1.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/'),
                           'energy_plot2': '{}{}energy2.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/'),
                           'energy_plot3': '{}{}energy3.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/')
                           }

        template_calc_specific = get_calc_specific_latex_template(
            TexTempl.CylindricalRollerBearing.value, sim)
        with open(calc_spec_tex_file_handle, 'w') as f:
            f.write(template_calc_specific.render(latex_variables))

    def generate_latex_figures(self, ui=None, res_dir=None):
        """Append system-specific output data to latex file"""
        max_press_field_ring1 = self.ring1.press[:, :, 0, 0]
        max_press_field_ring2 = self.ring2.press[:, :, 0, 0]

        plt_profile_approx(res_dir, SubDir.tex_figs.value)
        plt_contact(self.roller, self.ring1, PltOpts.DDD.value, res_dir,
                    SubDir.tex_figs.value, 'contact1')
        plt_contact(self.roller, self.ring2, PltOpts.DDD.value, res_dir,
                    SubDir.tex_figs.value, 'contact2')
        polplt_scatt_line(self.pol_coords_slip, self.slip_data, self.pol_ax,
                          self.roller.slip, 'roller slip',
                          'roller slip', res_dir, SubDir.tex_figs.value,
                          'slip1')

        plt_3d(self.roller.x_axis, self.roller.y_axis, max_press_field_ring1,
               self.roller.x_label, self.roller.y_label,
               'Pressure in MPa', 'pressure at highest normal force', res_dir,
               SubDir.tex_figs.value,
               'pressure1')

        plt_3d(self.roller.x_axis, self.roller.y_axis, max_press_field_ring2,
               self.roller.x_label, self.roller.y_label,
               'Pressure in MPa', 'pressure at highest normal force', res_dir,
               SubDir.tex_figs.value,
               'pressure2')

        polplt_scatt_line(np.append(self.phi_mat[0, :], self.phi_mat[0, 0]),
                          np.append(self.roller_norm_forces[0, :],
                                    self.roller_norm_forces[0, 0]),
                          np.append(self.phi_mat[0, :], self.phi_mat[0, 0]),
                          np.append(self.roller_norm_forces[0, :],
                                    self.roller_norm_forces[0, 0]),
                          'normal force in N', 'load distribution', res_dir,
                          SubDir.tex_figs.value,
                          'load1')

        plot_axis_pol = np.linspace(-math.pi,
                                    math.pi - 2 * math.pi / self.res_pol,
                                    int(self.res_pol))
        plot_data_ring1 = np.concatenate(
            (self.ring1.pv[round(self.res_pol / 2) - 1:-1, :],
             self.ring2.pv[0:round(self.res_pol / 2), :]), axis=0)
        plot_data_ring2 = np.concatenate(
            (self.ring2.pv[round(self.res_pol / 2) - 1:-1, :],
             self.ring2.pv[0:round(self.res_pol / 2), :]), axis=0)

        polplt_line_line(plot_axis_pol,
                         plot_data_ring1[:, round(self.roller.res_x / 2) - 1],
                         plot_axis_pol,
                         plot_data_ring2[:, round(self.roller.res_x / 2) - 1],
                         'pv in {}'.format(Unit.pvrel.value),
                         'pv_rel inner (max) and outer ring',
                         'inner ring (max)', 'outer ring',
                         res_dir, SubDir.tex_figs.value, 'energy1')

        if self.rot_vel1 == 0:
            plot_data = np.concatenate((self.ring1.e_akin.transpose()[
                                        round(self.res_pol / 2) - 1:-1, :],
                                        self.ring1.e_akin.transpose()[
                                        0:round(self.res_pol / 2), :]), axis=0)
            polplt_line(plot_axis_pol,
                        plot_data[:, round(self.roller.res_x / 2) - 1],
                        'e_a,kin_rel max in {}'.format(Unit.eakin.value),
                        'e_a,kin max inner ring',
                        res_dir, SubDir.tex_figs.value, 'energy2')
        else:
            plt_2d_scatt_line(self.roller.x_axis, self.ring1.e_akin,
                              self.roller.x_axis, self.ring1.e_akin,
                              self.roller.x_label,
                              'average e_a,kin in {}'.format(Unit.eakin.value),
                              'average e_a,kin inner ring', res_dir,
                              SubDir.tex_figs.value, 'energy2')
        if self.rot_vel2 == 0:
            plot_data = np.concatenate(
                (self.ring2.e_akin.transpose()[round(self.res_pol / 2) - 1:-1,
                 :],
                 self.ring2.e_akin.transpose()[0:round(self.res_pol / 2), :]),
                axis=0)
            polplt_line(plot_axis_pol,
                        plot_data[:, round(self.roller.res_x / 2) - 1],
                        'e_a,kin_rel max in {}'.format(Unit.eakin.value),
                        'e_a,kin max outer ring',
                        res_dir, SubDir.tex_figs.value, 'energy3')
        else:
            plt_2d_scatt_line(self.roller.x_axis, self.ring2.e_akin,
                              self.roller.x_axis, self.ring2.e_akin,
                              self.roller.x_label,
                              'average e_a,kin in {}'.format(Unit.eakin.value),
                              'average e_a,kin outer ring', res_dir,
                              SubDir.tex_figs.value, 'energy3')

    def plt_profiles_contacts(self, res_dir):
        """Generate plots of profiles and contacts"""
        print_it("plotting contacts", PrintOpts.lvl1.value)
        plt_contact(self.roller, self.ring1, PltOpts.DD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.roller, self.ring1, PltOpts.DDD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.roller, self.ring2, PltOpts.DD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.roller, self.ring2, PltOpts.DDD.value, res_dir,
                    SubDir.contacts.value)

    def plt_load_distr(self, res_dir):
        """Generate plots of load distribution"""
        print_it("plotting load distributions", PrintOpts.lvl1.value)
        for roller_position in range(self.num_roller_pos):
            roller_forces = self.roller_norm_forces[roller_position, :]
            roller_forces_for_plot = np.append(roller_forces, roller_forces[0])
            phi_plot = np.append(self.phi_mat[roller_position, :],
                                 self.phi_mat[roller_position, 0])
            polplt_scatt_line(phi_plot, roller_forces_for_plot, phi_plot,
                              roller_forces_for_plot, 'normal force in N',
                              ('load distribution position ' + str(
                                  roller_position)), res_dir,
                              SubDir.load_distr.value, (
                              'load_distribution_position_' + str(
                                  roller_position)),
                              self.roller_norm_forces[0, 0] * 1.1)
            print_progress(roller_position, self.num_roller_pos)

        polplt_scatt_line(self.pol_ax, self.normal_forces, self.pol_ax,
                          self.normal_forces, 'normal force in N',
                          'roller normal force', res_dir,
                          SubDir.load_distr.value, 'normal-force-over-position')

    def plt_pol_slip(self, res_dir):
        """Generate plot of polar roller slip"""
        polplt_scatt_line(self.pol_coords_slip, self.slip_data, self.pol_ax,
                          self.roller.slip, 'roller slip',
                          'roller slip', res_dir, SubDir.energy.value,
                          'roller.slip')

    def plt_pres(self, res_dir):
        """Generate contact pressure plots"""
        print_it("plotting pressure distributions", PrintOpts.lvl1.value)
        press_count = 0
        z_limit = np.amax(
            self.ring1.uniq_press[:, :, self.uniq_norm_forces.shape[0] - 1])
        for uniq_norm_force in range(self.uniq_norm_forces.shape[0] - 1, -1,
                                     -1):
            plt_3d(self.roller.x_axis, self.roller.y_axis,
                   self.ring1.uniq_press[:, :, uniq_norm_force],
                   self.roller.x_label, self.roller.y_label, 'pressure in MPa',
                   ('contact_pressure_ir_' + str(int(press_count))), res_dir,
                   SubDir.pressures.value,
                   ('contact_pressure_ir_' + str(int(press_count))),
                   [0, z_limit])
            plt_3d(self.roller.x_axis, self.roller.y_axis,
                   self.ring2.uniq_press[:, :, uniq_norm_force],
                   self.roller.x_label, self.roller.y_label, 'pressure in MPa',
                   ('contact_pressure_or_' + str(int(press_count))), res_dir,
                   SubDir.pressures.value,
                   ('contact_pressure_or_' + str(int(press_count))),
                   [0, z_limit])
            press_count = print_progress(press_count,
                                         self.uniq_norm_forces.shape[0])

    def plt_pv(self, res_dir):
        """Plot pv"""
        plt_ax_pol = np.linspace(-math.pi, math.pi - 2 * math.pi / self.res_pol,
                                 int(self.res_pol))
        plot_data_ring1 = np.concatenate(
            (self.ring1.pv[round(self.res_pol / 2) - 1:-1, :],
             self.ring2.pv[0:round(self.res_pol / 2), :]), axis=0)
        plot_data_ring2 = np.concatenate(
            (self.ring2.pv[round(self.res_pol / 2) - 1:-1, :],
             self.ring2.pv[0:round(self.res_pol / 2), :]), axis=0)
        plt_3d(plt_ax_pol, self.roller.x_axis, plot_data_ring1,
               'polar ring coordinate', 'roller length in mm',
               'pv in {}'.format(Unit.pvrel.value), 'pv inner ring', res_dir,
               SubDir.energy.value, 'max_pv_inner_ring',
               None, PltOpts.azim.value)
        plt_3d(plt_ax_pol, self.roller.x_axis, plot_data_ring2,
               'polar ring coordinate', 'roller length in mm',
               'pv in {}'.format(Unit.pvrel.value), 'pv outer ring', res_dir,
               SubDir.energy.value, 'pv_outer_ring',
               None, PltOpts.azim.value)

        polplt_line_line(plt_ax_pol,
                         plot_data_ring1[:, round(self.roller.res_x / 2) - 1],
                         plt_ax_pol,
                         plot_data_ring2[:, round(self.roller.res_x / 2) - 1],
                         'pv in {}'.format(Unit.pvrel.value),
                         'pv_rel inner (max) and outer ring',
                         'inner ring (max)', 'outer ring', res_dir,
                         SubDir.energy.value, 'pv-rel-max')
        scipy.io.savemat(res_dir +
                         '{}{}{}pv_rel_ring2.mat'.format(
                             os.sep,
                             SubDir.matlab.value,
                             os.sep
                         ),
                         dict(plot_axis_pol=plt_ax_pol,
                              plot_axis_x=self.roller.x_axis,
                              plot_data_pv_rel_ir_3d=plot_data_ring1,
                              plot_data_pv_rel_ir_2d=np.amax(plot_data_ring1,
                                                             axis=1),
                              plot_data_pv_rel_or_3d=plot_data_ring2,
                              plot_data_pv_rel_or_2d=np.amax(plot_data_ring2,
                                                             axis=1),
                              axis_pol=self.pol_ax,
                              axis_x=self.roller.x_axis,
                              data_pv_rel_or=self.ring2.pv))

    def plt_eakin(self, res_dir):
        """Plot kinetic friction energy accumulation"""
        print_it("plotting energy figures", PrintOpts.lvl1.value)
        # e_a,kin
        plt_ax_pol = np.linspace(-math.pi, math.pi - 2 * math.pi / self.res_pol,
                                 int(self.res_pol))
        if self.rot_vel1 == 0:
            plot_data = np.concatenate((self.ring1.e_akin.transpose()[
                                        round(self.res_pol / 2) - 1:-1, :],
                                        self.ring1.e_akin.transpose()[
                                        0:round(self.res_pol / 2), :]), axis=0)
            plt_3d(plt_ax_pol, self.roller.x_axis, plot_data,
                   'polar ring coordinate', 'roller length in mm',
                   'e_a,kin in {}'.format(Unit.eakin.value),
                   'e_a,kin inner ring', res_dir, SubDir.energy.value,
                   'e-a-kin-inner-ring', None, PltOpts.azim.value)

            polplt_line(plt_ax_pol,
                        plot_data[:, round(self.roller.res_x / 2) - 1],
                        'e_a,kin_rel max in {}'.format(Unit.eakin.value),
                        'e_a,kin max inner ring', res_dir,
                        SubDir.energy.value, 'max-e-a-kin-inner-ring')
        else:
            plt_2d_scatt_line(self.roller.x_axis, self.ring1.e_akin,
                              self.roller.x_axis, self.ring1.e_akin,
                              self.roller.x_label,
                              'average e_a,kin in {}'.format(Unit.eakin.value),
                              'average e_a,kin inner ring', res_dir,
                              SubDir.energy.value, 'average-e-a-kin-inner-ring')

        if self.rot_vel2 == 0:
            plot_data = np.concatenate((self.ring2.e_akin.transpose()[
                                        round(self.res_pol / 2) - 1:-1, :],
                                        self.ring2.e_akin.transpose()[
                                        0:round(self.res_pol / 2), :]), axis=0)
            plt_3d(plt_ax_pol, self.roller.x_axis, plot_data,
                   'polar ring coordinate', 'roller length in mm',
                   'e_a,kin in {}'.format(Unit.eakin.value),
                   'e_a,kin outer ring', res_dir, SubDir.energy.value,
                   'e-a-kin-outer-ring', None, PltOpts.azim.value)

            polplt_line(plt_ax_pol,
                        plot_data[:, round(self.roller.res_x / 2) - 1],
                        'e_a,kin_rel max in {}'.format(Unit.eakin.value),
                        'e_a,kin max outer ring', res_dir,
                        SubDir.energy.value, 'max-e-a-kin-outer-ring')
        else:
            plt_2d_scatt_line(self.roller.x_axis, self.ring2.e_akin,
                              self.roller.x_axis, self.ring2.e_akin,
                              self.roller.x_label,
                              'average e_a,kin in {}'.format(Unit.eakin.value),
                              'average e_a,kin outer ring', res_dir,
                              SubDir.energy.value, 'average-e-a-kin-outer-ring')

    def generate_summary_figure(self, res_dir):
        """Generate figure summarising the simulation results in a polar plot"""
        plot_axis_pol = np.linspace(0, 2 * math.pi - 2 * math.pi / self.res_pol,
                                    int(self.res_pol))
        pv_max_ring1 = self.ring1.pv
        eakin_dat = self.ring1.e_akin.transpose()
        press_dat = self.ring1.max_press
        force_dat = self.normal_forces
        slip_dat = self.roller.slip

        plt.rc('grid', linewidth=0.5, linestyle=':')
        plt.figure(facecolor='w')
        ax = plt.subplot(111, polar=True)

        ax.grid(True)

        # plot pv
        r = np.amax(pv_max_ring1, axis=1) / np.amax(pv_max_ring1)
        theta = plot_axis_pol
        leg = 'pv'
        ax.plot(theta, r, 'k-', label=leg)
        ind = round(len(theta) * 0.17)
        thisr, thistheta = r[ind], theta[ind]
        ax.annotate(leg,
                    xy=(thistheta, thisr),  # theta, radius
                    xytext=(thistheta, thisr),  # fraction, fraction
                    textcoords='data',
                    arrowprops=dict(facecolor='black', shrink=0.05),
                    horizontalalignment='left',
                    verticalalignment='bottom',
                    )

        # plot e_akin
        r = np.amax(eakin_dat, axis=1) / np.amax(eakin_dat)
        leg = '$\mathregular{e_{a,kin}}$'
        ax.plot(theta, r, 'k--', label=leg)
        ind = round(len(theta) * 0.03)
        thisr, thistheta = r[ind], theta[ind]
        ax.annotate(leg,
                    xy=(thistheta, thisr),  # theta, radius
                    xytext=(thistheta, thisr),  # fraction, fraction
                    textcoords='data',
                    arrowprops=dict(facecolor='black', shrink=0.05),
                    horizontalalignment='right',
                    verticalalignment='top',
                    )

        # plot pressure
        r = np.amax(press_dat, axis=1) / np.amax(press_dat)
        leg = '$\mathregular{p_{max}}$'
        ax.plot(theta, r, 'k-.', label=leg)
        ind = round(len(theta) / 1.1)
        thisr, thistheta = r[ind], theta[ind]
        ax.annotate(leg,
                    xy=(thistheta, thisr),  # theta, radius
                    xytext=(thistheta, thisr),  # fraction, fraction
                    textcoords='data',
                    arrowprops=dict(facecolor='black', shrink=0.05),
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    )

        # plot force
        r = np.amax(force_dat, axis=1) / np.amax(force_dat)
        leg = '$\mathregular{F_{N}}$'
        ax.plot(theta, r, 'k--', label=leg)
        ind = round(len(theta) / 1.05)
        thisr, thistheta = r[ind], theta[ind]
        ax.annotate(leg,
                    xy=(thistheta, thisr),  # theta, radius
                    xytext=(thistheta, thisr),  # fraction, fraction
                    textcoords='data',
                    arrowprops=dict(facecolor='black', shrink=0.05),
                    horizontalalignment='left',
                    verticalalignment='bottom',
                    )

        # plot slip
        r = slip_dat / np.amax(slip_dat)
        leg = 'roller slip'
        ax.plot(theta, r, 'k--', label=leg)
        ind = round(len(theta) * 0.79)
        thisr, thistheta = r[ind], theta[ind]
        ax.annotate(leg,
                    xy=(thistheta, thisr),  # theta, radius
                    xytext=(thistheta, thisr),  # fraction, fraction
                    textcoords='data',
                    arrowprops=dict(facecolor='black', shrink=0.05),
                    horizontalalignment='left',
                    verticalalignment='bottom',
                    )

        # plt.legend(loc='upper center')
        ax.set_theta_zero_location("S")
        ax.set_rlabel_position(112.5)
        label_position = ax.get_rlabel_position()
        r_max = 1.1
        r_label = 'normalised value'
        ax.text(np.radians(label_position + 15), r_max / 2, r_label,
                rotation=22.5, ha='center', va='center')
        ax.set_rmax(r_max)
        make_directory(res_dir, SubDir.pub_figs.value)
        file_handle = res_dir + os.sep + SubDir.pub_figs.value + os.sep + \
                      'normalised-inner-ring' + '.png'
        plt.savefig(file_handle, bbox_inches='tight', dpi=300)
        plt.close()
