import math
import os
import re

import numpy as np

from BluPrintTriboSys import TriboSys
from Constants import PltOpts, SubDir, TexTempl, UnitTex, Unit, PrintOpts, \
    PreSol
from cartesian_plot_functions import plt_profile, plt_contact, plt_3d, \
    plt_2d_scatt_line, \
    plt_energy_ring_on_ring, plt_profile_approx
from generate_latex_output import get_calc_specific_latex_template
from hertz_equations import hertz_displ
from influ_matrix_management import load_influ_mat, cache_influ_mat
from solve_half_space import solve_half_space, pre_solve_half_space
from system_functions import print_it, to_preci, exit_program, save_to_matlab


class RingOnRing(TriboSys):
    """Global tribosystem data"""

    def __init__(self, num_planets, sun, planet, global_force,
                 setup_name='ring-on-ring'):
        super().__init__(setup_name, global_force, None)

        self.num_planets = num_planets
        self.sun = sun
        self.planet = planet
        self.sun.norm_forces = None
        self.init_force = None

        self.planet_slip = None
        self.rot_velocity = None
        self.rot_velocity2 = None
        self.sun_rot_vel = None
        self.planet_rot_vel = None
        self.sliding_vel = None
        self.sun.press = None
        self.sun.max_press = None

        self.sun.rot_vel = None
        self.planet.rot_vel = None
        self.sun.omega = None
        self.planet.omega = None
        self.sun.vel = None
        self.planet.vel = None
        self.slip = None
        self.rel_vel = None
        self.pv = None

        self.influ_mat_db_1 = None
        self.press_zone_len = None

    def calc_load_distribution(self, ui=None, res_dir=None):
        """Calculate load distribution"""
        print_it("calculating load distribution")
        self.sun.norm_forces = np.multiply(np.ones(self.num_planets),
                                           self.global_force / self.num_planets)
        self.init_force = self.global_force

    def get_grid_size(self, ui, res_dir):
        """Determine grid size by running (quick) simulation with simplified
        contact bodies"""
        print_it('determining grid size', PrintOpts.lvl1.value)
        self.sun.simple_clone()
        self.sun.clone.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                                    self.init_force)
        self.planet.simple_clone()
        self.planet.clone.make_slave_to(self.sun.clone)

        init_displ = hertz_displ(self.sun.clone.e, self.planet.e,
                                 self.sun.clone.ny, self.planet.ny,
                                 self.sun.clone.r_hertz_x,
                                 self.sun.clone.r_hertz_y,
                                 self.planet.clone.r_hertz_x,
                                 self.planet.clone.r_hertz_y,
                                 self.sun.norm_forces[0])
        too_many_els_in_y = 1
        too_many_els_in_x = 1
        contact_width_y = 0.05
        contact_width_x = 0.05
        while too_many_els_in_y != 0 or \
              too_many_els_in_x != 0:
            self.sun.clone.make_profile(self.sun.clone.res_x,
                                        self.sun.clone.res_y, self.init_force,
                                        contact_width=contact_width_y,
                                        contact_length=contact_width_x)
            self.planet.clone.make_slave_to(self.sun.clone)

            pressure, init_displ = \
                pre_solve_half_space(self.sun.clone.profile,
                                     self.planet.clone.profile,
                                     self.sun.clone.x_axis,
                                     self.sun.clone.y_axis,
                                     self.sun.clone.res_x, self.sun.clone.res_y,
                                     self.sun.clone.delta_x,
                                     self.sun.clone.delta_y, self.sun.clone.e,
                                     self.planet.clone.e, self.sun.clone.ny,
                                     self.planet.clone.ny,
                                     self.sun.norm_forces[0],
                                     init_displ=init_displ, print_prog=False)

            pressure_els_y = sum(
                pressure[math.floor(self.sun.clone.res_y / 2), :] > 0)
            too_many_els_in_y = self.sun.clone.res_y - pressure_els_y - 2
            if too_many_els_in_y:
                contact_width_y += -np.sign(
                    too_many_els_in_y) * contact_width_y / 25

            pressure_els_x = sum(
                pressure[:, math.floor(self.sun.clone.res_x / 2)] > 0)
            too_many_els_in_x = self.sun.clone.res_x - pressure_els_x - 2
            if too_many_els_in_x:
                contact_width_x += -np.sign(
                    too_many_els_in_x) * contact_width_x / 25

        self.sun.make_profile(self.sun.res_x, self.sun.res_y, self.init_force,
                              contact_width=contact_width_y,
                              contact_length=contact_width_x)
        self.planet.make_slave_to(self.sun)
        return init_displ

    def calc_contact_pressure(self, ui=None, res_dir=None):
        """Calculate contact pressure distribution between sun and planet
        ring(s)"""
        print_it('calculating 1 pressure distribution')
        init_displ = self.get_grid_size(ui, res_dir)
        [self.influ_mat_db_1] = load_influ_mat(ui, res_dir, 1)
        print_it('solving first half space', PrintOpts.lvl1.value)
        self.sun.press, self.influ_mat_db_1 = \
            solve_half_space(self.sun.profile, self.planet.profile,
                             self.sun.x_axis, self.sun.y_axis, self.sun.res_x,
                             self.sun.res_y, self.sun.delta_x, self.sun.delta_y,
                             self.sun.e, self.planet.e,
                             self.sun.ny, self.planet.ny,
                             self.sun.norm_forces[0], res_dir,
                             init_displ=init_displ,
                             influ_mat_db=self.influ_mat_db_1)

        cache_influ_mat(ui, [self.influ_mat_db_1], res_dir)
        self.sun.max_press = np.amax(self.sun.press, axis=1)

        dat_dict = dict(x_axis=self.sun.x_axis,
                        y_axis=self.sun.y_axis,
                        contact_pressure=self.sun.press)
        save_to_matlab(dat_dict, res_dir, 'pressure_field')

    def calc_kinematics(self, rot_vel1, rot_vel2, ui=None, res_dir=None):
        """Calculate tribosystem kinematics based on rotational velocities of
        sun and planet(s)"""
        print_it("calculating kinematics")
        self.sun.rot_vel = rot_vel1
        self.planet.rot_vel = rot_vel2
        self.sun.omega = self.sun.rot_vel / 60
        self.planet.omega = self.planet.rot_vel / 60
        self.sun.vel = self.sun.diameter * math.pi * self.sun.omega
        self.planet.vel = self.planet.diameter * math.pi * self.planet.omega
        self.slip = (self.sun.vel - self.planet.vel) / self.sun.vel
        self.rel_vel = np.ones(self.sun.res_x) * (
        self.sun.vel - self.planet.vel)

        self.sun.footpr_vel = \
            2 * math.pi * self.sun.diameter / 2 * self.sun.omega
        self.planet.footpr_vel = \
            2 * math.pi * self.planet.diameter / 2 * self.planet.omega
        try:
            self.sun.overroll_t_incr = self.sun.delta_y / self.sun.footpr_vel
            self.planet.overroll_t_incr = \
                self.planet.delta_y / self.planet.footpr_vel
        except ZeroDivisionError:
            exit_program(
                'rotational velocities of sun and planet must not be 0')

        self.press_zone_len = (self.sun.press > 0).sum(1) * self.sun.delta_y
        self.sun.overroll_t = np.divide(self.press_zone_len,
                                        self.sun.footpr_vel)
        self.planet.overroll_t = np.divide(self.press_zone_len,
                                           self.planet.footpr_vel)
        self.sun.no_overroll_t = np.divide(
            (2 * math.pi * (self.sun.diameter / 2) - self.num_planets *
             self.press_zone_len) / self.num_planets, self.sun.footpr_vel)
        self.planet.no_overroll_t = np.divide(
            (2 * math.pi * (self.planet.diameter / 2) - self.press_zone_len),
            self.planet.footpr_vel)

    def calc_pv(self, ui=None, res_dir=None):
        """Calculate product of local maximum pressure and local maximum
        relative velocity"""
        print_it("calculating pv_rel")
        self.pv = np.multiply(abs(self.rel_vel), self.sun.max_press) / 1000
        dat_dict = dict(x_axis=self.sun.x_axis,
                        pv_rel=self.pv)
        save_to_matlab(dat_dict, res_dir, 'pv-rel')

    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation in W per m^2"""
        print_it("calculating e_a,kin")
        pv_local = np.multiply(self.sun.press.sum(1), self.rel_vel)
        self.sun.e_akin = np.absolute(
            np.divide(np.multiply(pv_local, self.sun.overroll_t_incr),
                      self.sun.no_overroll_t)) / 1000
        self.planet.e_akin = np.absolute(
            np.divide(np.multiply(pv_local, self.planet.overroll_t_incr),
                      self.planet.no_overroll_t)) / 1000

    def plot_it(self, ui=None, res_dir=None):
        """Orchestrate output plot generation"""
        print_it("plotting results")
        plt_profile(self.sun, PltOpts.DD.value, res_dir, SubDir.profiles.value)
        plt_profile(self.sun, PltOpts.DDD.value, res_dir, SubDir.profiles.value)
        plt_profile(self.planet, PltOpts.DD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.planet, PltOpts.DDD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile_approx(res_dir, SubDir.profiles.value)
        plt_contact(self.sun, self.planet, PltOpts.DD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.sun, self.planet, PltOpts.DDD.value, res_dir,
                    SubDir.contacts.value)

        plt_3d(self.sun.x_axis, self.sun.y_axis, self.sun.press,
               self.sun.x_label, self.sun.y_label, 'pressure in MPa',
               'contact_pressure_sun', res_dir, SubDir.pressures.value,
               'contact_pressure_sun')
        plt_2d_scatt_line(self.sun.x_axis, self.pv, self.sun.x_axis, self.pv,
                          self.sun.x_label,
                          'pv_rel in {}'.format(Unit.pvrel.value), 'pv_rel',
                          res_dir, SubDir.energy.value, 'pv_rel')
        plt_2d_scatt_line(self.sun.x_axis, self.sun.e_akin, self.sun.x_axis,
                          self.sun.e_akin, self.sun.x_label,
                          'e_akin in {}'.format(Unit.eakin.value), 'e_akin',
                          res_dir, SubDir.energy.value, 'sun.e_akin')
        plt_2d_scatt_line(self.planet.x_axis, self.planet.e_akin,
                          self.planet.x_axis, self.planet.e_akin,
                          self.planet.x_label,
                          'e_akin in {}'.format(Unit.eakin.value), 'e_akin',
                          res_dir,
                          SubDir.energy.value, 'planet.e_akin')
        plt_energy_ring_on_ring(self, res_dir, SubDir.energy.value,
                                'e-akin-vs-pv-rel')

    def generate_latex_output(self, calc_spec_tex_file_handle, sim, ui=None,
                              res_dir=None):
        """Generate calculation-specific part of the LaTeX output file"""
        average_pressure = np.mean(self.sun.press[self.sun.press > 0])
        numeric_output_data = [
            ('pressure, max.', to_preci(np.amax(self.sun.press), 4),
             UnitTex.pressure.value, 'unverified'),
            ('pressure, av.', to_preci(average_pressure, 4),
             UnitTex.pressure.value, 'unverified'),
            ('e_a,kin sun, max.', to_preci(np.amax(self.sun.e_akin), 4),
             UnitTex.eakin.value, 'unverified'),
            ('e_a,kin planet, max.', to_preci(np.amax(self.planet.e_akin), 4),
             UnitTex.eakin.value, 'unverified'),
            ('pv_rel, max.', to_preci(np.amax(self.pv), 4),
             UnitTex.pvrel.value, 'unverified'),
            ('contact area', to_preci(self.sun.get_area(self.sun.press), 4),
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
            TexTempl.RingOnRing.value, sim)
        with open(calc_spec_tex_file_handle, 'w') as f:
            f.write(template_calc_specific.render(latex_variables))

    def generate_latex_figures(self, ui=None, res_dir=None):
        """Generate calculation-specific figures for LaTeX report"""
        plt_contact(self.sun, self.planet, PltOpts.DDD.value, res_dir,
                    SubDir.tex_figs.value, 'contact1')
        plt_profile_approx(res_dir, SubDir.tex_figs.value)
        plt_3d(self.sun.x_axis, self.sun.y_axis, self.sun.press,
               self.sun.x_label, self.sun.y_label, 'pressure in MPa',
               'contact_pressure_sun', res_dir, SubDir.tex_figs.value,
               'pressure1')
        plt_energy_ring_on_ring(self, res_dir, SubDir.tex_figs.value, 'energy1')
