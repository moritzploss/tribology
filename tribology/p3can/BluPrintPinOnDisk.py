import math
import os
import re

import numpy as np

from BluPrintTriboSys import TriboSys
from Constants import PltOpts, SubDir, TexTempl, Unit, UnitTex, PrintOpts, \
    PreSol
from cartesian_plot_functions import plt_contact, plt_3d, plt_2d, \
    plt_2d_scatt_line, plt_2d_2y_ax, \
    plt_profile_approx, plt_cum_dist, plt_profiles
from generate_latex_output import get_calc_specific_latex_template
from hertz_equations import hertz_displ
from influ_matrix_management import load_influ_mat, cache_influ_mat
from solve_half_space import solve_half_space, pre_solve_half_space
from system_functions import print_it, to_preci, save_to_matlab


class PinOnDisk(TriboSys):
    # TODO: add class attribute description
    def __init__(self, num_pins, sliding_diameter, pin, disk, global_force,
                 setup_name='pin-on-disk'):
        super().__init__(setup_name, global_force, None)

        self.num_pins = num_pins
        self.sliding_diameter = sliding_diameter  # in mm
        self.raceway_vel = None
        self.press_zone_len = None
        self.effective_rot_velocity = None
        self.effective_omega = None

        self.influ_mat_db_1 = None

        self.rot_vel1 = None
        self.norm_forces = None

        self.pin = pin
        self.pin.press = None
        self.pin.max_press = None
        self.pin.rel_vel = None

        self.disk = disk
        self.disk.press = None
        self.disk.max_press = None

        self.pv = None
        self.e_akin = None

        self.footpr_vel = None
        self.overroll_t_incr = None
        self.overroll_t = None
        self.no_overroll_t = None

        if self.pin.diameter == self.pin.circle_radius:
            self.pin_is_ball = True
        else:
            self.pin_is_ball = False

    def calc_load_distribution(self, ui=None, res_dir=None):
        """Calculate load distribution"""
        print_it("calculating load distribution")
        self.norm_forces = np.multiply(np.ones(self.num_pins),
                                       self.global_force / self.num_pins)
        self.init_force = self.global_force / self.num_pins

    def get_grid_size(self, ui=None, res_dir=None):
        """Determine grid size by running (quick) simulation with simplified
        contact bodies"""
        print_it('determining grid size', PrintOpts.lvl1.value)
        self.pin.simple_clone()
        self.pin.clone.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                                    self.init_force)
        self.disk.simple_clone()
        self.disk.clone.make_slave_to(self.pin.clone)

        init_displ = hertz_displ(self.pin.clone.e, self.disk.clone.e,
                                 self.pin.clone.ny, self.disk.clone.ny,
                                 self.pin.clone.r_hertz_x,
                                 self.pin.clone.r_hertz_y,
                                 self.disk.clone.r_hertz_x,
                                 self.disk.clone.r_hertz_y, self.norm_forces[0])
        contact_delta = 1
        contact_width = 0.025
        while contact_delta != 0:
            if type(self.pin).__name__ == "Ring":
                self.pin.clone.make_profile(self.pin.clone.res_x,
                                            self.pin.clone.res_y,
                                            self.init_force,
                                            contact_width=contact_width)
            elif type(self.pin).__name__ == "Ball":
                self.pin.clone.make_profile(PreSol.res_x.value,
                                            PreSol.res_y.value, self.init_force,
                                            contact_width=contact_width)
                self.disk.make_slave_to(self.pin)

            pressure, init_displ = \
                pre_solve_half_space(self.pin.clone.profile,
                                     self.disk.clone.profile,
                                     self.pin.clone.x_axis,
                                     self.pin.clone.y_axis,
                                     self.pin.clone.res_x, self.pin.clone.res_y,
                                     self.pin.clone.delta_x,
                                     self.pin.clone.delta_y, self.pin.clone.e,
                                     self.disk.clone.e, self.pin.clone.ny,
                                     self.disk.clone.ny, self.norm_forces[0],
                                     init_displ=init_displ)
            pressure_elements = sum(
                pressure[math.floor(self.pin.clone.res_y / 2), :] > 0)
            contact_delta = pressure_elements + 2 - self.pin.clone.res_y
            contact_width += np.sign(contact_delta) * contact_width / 25

        if type(self.pin).__name__ == "Ring":
            self.pin.make_profile(self.pin.res_x, self.pin.res_y,
                                  self.init_force, contact_width=contact_width)
        elif type(self.pin).__name__ == "Ball":
            self.pin.make_profile(self.pin.res_x, self.pin.res_y,
                                  self.init_force, contact_width=contact_width)
            self.disk.make_slave_to(self.pin)
        return init_displ

    def calc_contact_pressure(self, ui=None, res_dir=None):
        """Calculate contact pressure distribution between pin and disk"""
        print_it('calculating 1 pressure distribution')

        init_displ = self.get_grid_size(ui, res_dir)
        [self.influ_mat_db_1] = load_influ_mat(ui, res_dir, 1)
        print_it('solving first half space', PrintOpts.lvl1.value)
        self.pin.press, self.influ_mat_db_1 = \
            solve_half_space(self.pin.profile, self.disk.profile,
                             self.pin.x_axis, self.pin.y_axis, self.pin.res_x,
                             self.pin.res_y, self.pin.delta_x, self.pin.delta_y,
                             self.pin.e, self.disk.e,
                             self.pin.ny, self.disk.ny, self.norm_forces[0],
                             res_dir, init_displ=init_displ,
                             influ_mat_db=self.influ_mat_db_1)

        cache_influ_mat(ui, [self.influ_mat_db_1], res_dir)
        self.pin.max_press = np.amax(self.pin.press, axis=1)
        self.save_pressure_to_matlab(res_dir)

    def save_pressure_to_matlab(self, res_dir):
        dat_dict = dict(contact_pressure=self.pin.press,
                        x_axis=self.pin.x_axis,
                        y_axis=self.pin.y_axis)
        save_to_matlab(dat_dict, res_dir, 'pressure-field')

    def calc_kinematics(self, rot_vel1, rot_vel2, ui=None, res_dir=None):
        """Calculate bearing kinematics based on rotational velocity of disk"""
        print_it("calculating kinematics")
        self.rot_vel1 = rot_vel1
        self.effective_rot_velocity = self.rot_vel1
        self.effective_omega = self.effective_rot_velocity / 60
        self.pin.rel_vel = np.add(-self.pin.x_axis,
                                  self.sliding_diameter) * \
                                        math.pi * self.effective_omega

        self.press_zone_len = (self.pin.press > 0).sum(1) * self.pin.delta_y
        self.footpr_vel = self.pin.rel_vel
        self.overroll_t_incr = self.pin.delta_y / self.footpr_vel
        self.overroll_t = np.divide(self.press_zone_len, self.footpr_vel)
        self.no_overroll_t = np.divide((2 * math.pi * (
        self.sliding_diameter / 2) - self.num_pins * self.press_zone_len)
                                       / self.num_pins, self.footpr_vel)

    def calc_pv(self, ui=None, res_dir=None):
        """Calculate product of local maximum pressure and local maximum
        relative velocity"""
        print_it("calculating pv_rel")
        self.pv = np.multiply(self.pin.max_press, self.pin.rel_vel) / 1000

    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation in W per m^2"""
        print_it("calculating e_a,kin")
        pv_local = np.multiply(self.pin.press.sum(1), self.pin.rel_vel)
        self.e_akin = np.absolute(
            np.divide(np.multiply(pv_local, self.overroll_t_incr),
                      self.no_overroll_t)) / 1000

    def plot_it(self, ui=None, res_dir=None):
        """Orchestrate output plot generation"""
        print_it("plotting results")
        plt_profiles(res_dir)
        print_it("plotting other", PrintOpts.lvl1.value)
        plt_contact(self.pin, self.disk, PltOpts.DD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.pin, self.disk, PltOpts.DDD.value, res_dir,
                    SubDir.contacts.value)
        plt_3d(self.pin.x_axis, self.pin.y_axis, self.pin.press,
               self.pin.x_label, self.pin.y_label,
               'pressure in MPa', 'contact_pressure_pin', res_dir,
               SubDir.pressures.value, 'contact_pressure')
        plt_2d(self.pin.x_axis, self.pin.rel_vel, self.pin.x_label, 'rel vel',
               'rel vel', res_dir, SubDir.energy.value,
               'rel_vel')
        plt_2d(self.pin.x_axis, self.pv, self.pin.x_label,
               'pv_rel in {}'.format(Unit.pvrel.value), 'pv_rel', res_dir,
               SubDir.energy.value, 'pv_rel')
        plt_2d_scatt_line(self.pin.x_axis, self.e_akin, self.pin.x_axis,
                          self.e_akin, self.pin.x_label,
                          'e_akin in W $m^{-2}$', 'e_akin', res_dir,
                          SubDir.energy.value, 'e_akin')
        plt_2d_2y_ax(self.pin.x_axis, self.e_akin, self.pin.x_axis, self.pv,
                     self.pin.x_label, 'e_akin vs pv_rel',
                     'e_akin in {}'.format(Unit.eakin.value),
                     'pv_rel in {}'.format(Unit.eakin.value),
                     res_dir, SubDir.energy.value, 'e_akin_vs_pv_rel')

        plt_cum_dist(self.pin.press, 'pressure in {}'.format(Unit.press.value),
                     res_dir, SubDir.pressures.value,
                     'cumulative pressure distribution')

    def generate_latex_output(self, calc_spec_tex_file_handle, sim, ui=None,
                              res_dir=None):
        """Generate calculation-specific part of the LaTeX output file"""
        average_pressure = np.mean(self.pin.press[self.pin.press > 0])
        numeric_output_data = [
            ('pressure, max.', to_preci(np.amax(self.pin.press), 4),
             UnitTex.pressure.value, 'unverified'),
            ('pressure, av.', to_preci(average_pressure, 4),
             UnitTex.pressure.value, 'unverified'),
            ('e_a,kin, max.', to_preci(np.amax(self.e_akin), 4),
             UnitTex.eakin.value, 'unverified'),
            ('pv_rel, max.', to_preci(np.amax(self.pv), 4),
             UnitTex.pvrel.value, 'unverified'),
            ('contact area', to_preci(self.pin.get_area(self.pin.press), 4),
             UnitTex.area.value, 'unverified')]

        table_calc_summary = []
        for key, value, unit, status in sorted(numeric_output_data):
            table_calc_summary.append(
                (re.sub('_', '\_', key), value, unit, status))

        latex_variables = {'table_calc_summary': table_calc_summary,
                           'contact_plot1': '{}{}contact1.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/'),
                           'pressure_plot1': '{}{}pressure1.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/'),
                           'energy_plot1': '{}{}energy1.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value, '/')}
        template_calc_specific = get_calc_specific_latex_template(
            TexTempl.PinOnDisk.value, sim)
        with open(calc_spec_tex_file_handle, 'w') as f:
            f.write(template_calc_specific.render(latex_variables))

    def generate_latex_figures(self, ui=None, res_dir=None):
        """Generate calculation-specific figures for LaTeX report"""
        plt_contact(self.pin, self.disk, PltOpts.DDD.value, res_dir,
                    SubDir.tex_figs.value, 'contact1')
        plt_profile_approx(res_dir, SubDir.tex_figs.value)
        plt_3d(self.pin.x_axis, self.pin.y_axis, self.pin.press,
               self.pin.x_label, self.pin.y_label, 'pressure in MPa',
               'contact_pressure_pin', res_dir, SubDir.tex_figs.value,
               'pressure1')

        plt_2d_2y_ax(self.pin.x_axis, self.e_akin, self.pin.x_axis, self.pv,
                     self.pin.x_label, 'e_akin vs pv_rel',
                     'e_akin in W $\mathregular{m^{-2}}$',
                     'pv_rel in W $\mathregular{m^{-2}}$', res_dir,
                     SubDir.tex_figs.value, 'energy1')
