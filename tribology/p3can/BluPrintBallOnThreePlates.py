import math
import os
import re

import numpy as np

from BluPrintTriboSys import TriboSys
from Constants import PltOpts, SubDir, TexTempl, UnitTex, Unit, PrintOpts, \
    PreSol
from cartesian_plot_functions import plt_profile, \
    plt_contact, plt_3d, plt_2d, plt_2d_scatt_line, plt_2d_2y_ax, \
    plt_profile_approx
from generate_latex_output import get_calc_specific_latex_template
from hertz_equations import hertz_displ
from influ_matrix_management import cache_influ_mat, load_influ_mat
from solve_half_space import solve_half_space, pre_solve_half_space
from system_functions import print_it, to_preci, save_to_matlab


class BallOnThreePlates(TriboSys):
    """Global tribosystem data"""

    def __init__(self, ball, plate, global_force, plate_angle=45,
                 setup_name='ball-on-3-plates'):
        super().__init__(setup_name, global_force, None)

        self.ball = ball
        self.plate = plate
        self.plate_angle = plate_angle / 180 * math.pi / 2

        self.sliding_diam = None
        self.rot_vel = None
        self.effective_omega = None
        self.raceway_vel = None
        self.ball.norm_forces = None
        self.ball.press = None
        self.ball.max_press = None
        self.rel_vel = None
        self.pv = None
        self.e_akin = None
        self.init_force = None

        self.influ_mat_db_1 = None
        self.footpr_vel = None
        self.overroll_t_incr = None
        self.press_zone_len = None
        self.overroll_t = None
        self.no_overroll_t = None

    def calc_load_distribution(self, ui=None, res_dir=None):
        """Calculate load distribution"""
        print_it("calculating load distribution")
        self.init_force = self.global_force / 3 / math.cos(self.plate_angle / 2)
        self.ball.norm_forces = np.ones(3) * self.init_force
        self.sliding_diam = 2 * self.ball.diameter / 2 * math.sin(
            (math.pi - self.plate_angle) / 2)

    def get_grid_size(self, ui=None, res_dir=None):
        """Determine grid size by running (quick) simulation with simplified
        contact bodies"""
        print_it('determining grid size', PrintOpts.lvl1.value)
        self.ball.simple_clone()
        self.ball.clone.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                                     self.init_force)
        self.plate.simple_clone()
        self.plate.clone.make_slave_to(self.ball.clone)

        init_displ = hertz_displ(self.ball.clone.e, self.plate.clone.e,
                                 self.ball.clone.ny, self.plate.clone.ny,
                                 self.ball.clone.r_hertz_x,
                                 self.ball.clone.r_hertz_y,
                                 self.plate.clone.r_hertz_x,
                                 self.plate.clone.r_hertz_y,
                                 self.ball.norm_forces[0])
        contact_delta = 1
        contact_width = 0.025
        while contact_delta != 0:
            self.ball.clone.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                                         self.init_force,
                                         contact_width=contact_width)
            self.plate.clone.make_slave_to(self.ball.clone)

            pressure, init_displ = \
                pre_solve_half_space(self.ball.clone.profile,
                                     self.plate.clone.profile,
                                     self.ball.clone.x_axis,
                                     self.ball.clone.y_axis,
                                     self.ball.clone.res_x,
                                     self.ball.clone.res_y,
                                     self.ball.clone.delta_x,
                                     self.ball.clone.delta_y, self.ball.clone.e,
                                     self.plate.clone.e, self.ball.clone.ny,
                                     self.plate.clone.ny,
                                     self.ball.norm_forces[0],
                                     init_displ=init_displ)

            pressure_elements = sum(
                pressure[math.floor(self.ball.clone.res_y / 2), :] > 0)
            contact_delta = pressure_elements + 2 - self.ball.clone.res_y
            contact_width += np.sign(contact_delta) * contact_width / 25

        self.ball.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                               self.init_force, contact_width=contact_width)
        self.plate.make_slave_to(self.ball)
        return init_displ

    def calc_contact_pressure(self, ui=None, res_dir=None):
        """Calculate contact pressure distribution between washer(s) and
        roller """
        print_it('calculating 1 pressure distribution')

        init_displ = self.get_grid_size(ui, res_dir)
        [self.influ_mat_db_1] = load_influ_mat(ui, res_dir, 1)
        print_it('solving first half space', PrintOpts.lvl1.value)
        self.ball.press, self.influ_mat_db_1 = \
            solve_half_space(self.ball.profile, self.plate.profile,
                             self.ball.x_axis, self.ball.y_axis,
                             self.ball.res_x,
                             self.ball.res_y, self.ball.delta_x,
                             self.ball.delta_y, self.ball.e,
                             self.plate.e,
                             self.ball.ny, self.plate.ny,
                             self.ball.norm_forces[0], res_dir,
                             init_displ=init_displ,
                             influ_mat_db=self.influ_mat_db_1)

        cache_influ_mat(ui, [self.influ_mat_db_1], res_dir)
        self.ball.max_press = np.amax(self.ball.press, axis=1)
        self.save_pressure_to_matlab(res_dir)

    def save_pressure_to_matlab(self, res_dir):
        dat_dict = dict(contact_pressure=self.ball.press,
                        x_axis=self.ball.x_axis,
                        y_axis=self.ball.y_axis)
        save_to_matlab(dat_dict, res_dir, 'pressure-field')

    def calc_kinematics(self, rot_vel1, rot_vel2, ui=None, res_dir=None):
        """Calculate tribosystem kinematics based on rotational velocities of
        contact bodies"""
        print_it("calculating pv_rel")
        self.effective_omega = rot_vel1 / 60
        self.raceway_vel = self.effective_omega * self.sliding_diam
        # neglecting sliding speed differences along ball contact area
        self.rel_vel = np.ones(self.ball.res_x) * self.raceway_vel

        self.footpr_vel = self.raceway_vel
        self.overroll_t_incr = np.divide(self.ball.delta_y, self.footpr_vel)
        self.press_zone_len = (self.ball.press > 0).sum(1) * self.ball.delta_y
        self.overroll_t = np.divide(self.press_zone_len, self.footpr_vel)
        self.no_overroll_t = np.divide((2 * math.pi * (
        self.sliding_diam / 2) - 3 * self.press_zone_len) / 3,
                                       self.footpr_vel)

    def calc_pv(self, ui=None, res_dir=None):
        """Calculate product of local maximum pressure and local maximum
        relative velocity"""
        print_it("calculating pv_rel")
        self.pv = np.multiply(self.ball.max_press, self.rel_vel) / 1000

    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation in W per m^2"""
        print_it("calculating e_a,kin")
        pv_local = np.multiply(self.ball.press.sum(1), self.rel_vel)
        self.e_akin = np.absolute(
            np.divide(np.multiply(pv_local, self.overroll_t_incr),
                      self.no_overroll_t)) / 1000

    def plot_it(self, ui=None, res_dir=None):
        """Orchestrate output plot generation"""
        print_it("plotting results and finishing up")
        plt_profile(self.ball, PltOpts.DD.value, res_dir, SubDir.profiles.value)
        plt_profile(self.ball, PltOpts.DDD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.plate, PltOpts.DD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.plate, PltOpts.DDD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile_approx(res_dir, SubDir.profiles.value)
        plt_contact(self.ball, self.plate, PltOpts.DD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.ball, self.plate, PltOpts.DDD.value, res_dir,
                    SubDir.contacts.value)

        plt_3d(self.ball.x_axis, self.ball.y_axis, self.ball.press,
               self.ball.x_label, self.ball.y_label,
               'pressure in MPa', 'contact pressure', res_dir,
               SubDir.pressures.value, 'contact-pressure')
        plt_2d(self.ball.x_axis, self.pv, self.ball.x_label,
               'pv_rel in {}'.format(Unit.pvrel.value), 'pv_rel',
               res_dir, SubDir.energy.value, 'pv_rel')
        plt_2d_scatt_line(self.ball.x_axis, self.e_akin, self.ball.x_axis,
                          self.e_akin, self.ball.x_label,
                          'e_akin in {}'.format(Unit.eakin.value), 'e_akin',
                          res_dir, SubDir.energy.value,
                          'e_akin')
        plt_2d_2y_ax(self.ball.x_axis, self.e_akin, self.ball.x_axis, self.pv,
                     self.ball.x_label, 'e_akin vs pv_rel',
                     'e_akin in {}'.format(Unit.eakin.value),
                     'pv_rel in {}'.format(Unit.pvrel.value),
                     res_dir, SubDir.energy.value, 'e_akin_vs_pv_rel')

    def generate_latex_output(self, calc_spec_tex_file_handle, sim, ui=None,
                              res_dir=None):
        """Generate the calculation-specific part of the latex output file"""
        average_pressure = np.mean(self.ball.press[self.ball.press > 0])
        numeric_output_data = [
            ('pressure, max.', to_preci(np.amax(self.ball.press), 4),
             UnitTex.pressure.value, 'unverified'),
            ('pressure, av.', to_preci(average_pressure, 4),
             UnitTex.pressure.value, 'unverified'),
            ('e_a,kin, max.', to_preci(np.amax(self.e_akin), 4),
             UnitTex.eakin.value, 'unverified'),
            ('pv_rel, max.', to_preci(np.amax(self.pv), 4),
             UnitTex.pvrel.value, 'unverified'),
            ('contact area', to_preci(self.ball.get_area(self.ball.press), 4),
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
            TexTempl.BallOn3Plates.value, sim)
        with open(calc_spec_tex_file_handle, 'w') as f:
            f.write(template_calc_specific.render(latex_variables))

    def generate_latex_figures(self, ui=None, res_dir=None):
        """Generate calculation-specific figures for LaTeX report"""
        plt_profile_approx(res_dir, SubDir.tex_figs.value)
        plt_contact(self.ball, self.plate, PltOpts.DDD.value, res_dir,
                    SubDir.tex_figs.value, 'contact1')
        plt_3d(self.ball.x_axis, self.ball.y_axis, self.ball.press,
               self.ball.x_label, self.ball.y_label,
               'pressure in MPa', 'contact pressure', res_dir,
               SubDir.tex_figs.value, 'pressure1')
        plt_2d_2y_ax(self.ball.x_axis, self.e_akin, self.ball.x_axis, self.pv,
                     self.ball.x_label, 'e_akin vs pv_rel',
                     'e_akin in {}'.format(UnitTex.eakin.value),
                     'pv_rel in {}'.format(UnitTex.pvrel.value),
                     res_dir, SubDir.tex_figs.value, 'energy1')
