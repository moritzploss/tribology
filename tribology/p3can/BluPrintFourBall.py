import math
import re

import numpy as np

from BluPrintTriboSys import TriboSys
from Constants import PltOpts, SubDir, TexTempl, UnitTex, PrintOpts, PreSol
from cartesian_plot_functions import plt_profile, plt_contact, plt_3d, \
    plt_2d_2y_ax, \
    plt_2d_scatt_line, plt_profile_approx
from generate_latex_output import get_calc_specific_latex_template
from hertz_equations import hertz_displ
from influ_matrix_management import cache_influ_mat, load_influ_mat
from solve_half_space import solve_half_space, pre_solve_half_space
from system_functions import print_it, to_preci, save_to_matlab


class FourBall(TriboSys):
    """Global tribosystem data"""

    def __init__(self, rot_ball, stat_ball, global_force,
                 setup_name='four-ball'):
        super().__init__(setup_name, global_force, None)

        self.sliding_diam = None  # in mm
        self.rot_ball = rot_ball
        self.stat_ball = stat_ball
        self.rot_vel = None

        self.effective_omega = None
        self.raceway_vel = None
        self.rot_ball.norm_forces = None
        self.rot_ball.max_press = None
        self.rel_vel = None
        self.pv = None
        self.e_akin = None
        self.init_force = None
        self.influ_mat_db_1 = None

        self.rot_ball.press = None
        self.rot_ball.max_press = None
        self.footpr_vel = None
        self.overroll_t_incr = None
        self.press_zone_len = None
        self.overroll_t = None
        self.no_overroll_t = None

    def calc_load_distribution(self, ui=None, res_dir=None):
        """Calculate load distribution"""
        print_it("calculating load distribution")
        r1 = self.rot_ball.diameter / 2
        r2 = self.stat_ball.diameter / 2
        r_circum_circle = math.sqrt(3) / 3 * 2 * r2
        contact_angle = math.acos(r_circum_circle / (r1 + r2))
        lever_arm = r_circum_circle - r2 * math.cos(contact_angle)
        self.sliding_diam = 2 * lever_arm
        self.init_force = self.global_force / math.sin(contact_angle) / 3
        self.rot_ball.norm_forces = np.ones(3) * self.init_force

    def get_grid_size(self, ui=None, res_dir=None):
        """Determine grid size by running (quick) simulation with simplified
        contact bodies"""
        print_it('determining grid size', PrintOpts.lvl1.value)
        self.rot_ball.simple_clone()
        self.rot_ball.clone.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                                         self.init_force)
        self.stat_ball.simple_clone()
        self.stat_ball.clone.make_slave_to(self.rot_ball.clone)

        init_displ = hertz_displ(self.rot_ball.clone.e, self.stat_ball.clone.e,
                                 self.rot_ball.clone.ny,
                                 self.stat_ball.clone.ny,
                                 self.rot_ball.clone.r_hertz_x,
                                 self.rot_ball.clone.r_hertz_y,
                                 self.stat_ball.clone.r_hertz_x,
                                 self.stat_ball.clone.r_hertz_y,
                                 self.rot_ball.norm_forces[0])
        contact_delta = 1
        contact_width = 0.025
        while contact_delta != 0:
            self.rot_ball.clone.make_profile(PreSol.res_x.value,
                                             PreSol.res_y.value,
                                             self.init_force,
                                             contact_width=contact_width)
            self.stat_ball.clone.make_slave_to(self.rot_ball.clone)

            pressure, init_displ = \
                pre_solve_half_space(self.rot_ball.clone.profile,
                                     self.stat_ball.clone.profile,
                                     self.rot_ball.clone.x_axis,
                                     self.rot_ball.clone.y_axis,
                                     self.rot_ball.clone.res_x,
                                     self.rot_ball.clone.res_y,
                                     self.rot_ball.clone.delta_x,
                                     self.rot_ball.clone.delta_y,
                                     self.rot_ball.clone.e,
                                     self.stat_ball.clone.e,
                                     self.rot_ball.clone.ny,
                                     self.stat_ball.clone.ny,
                                     self.rot_ball.norm_forces[0],
                                     init_displ=init_displ)

            pressure_elements = sum(
                pressure[math.floor(self.rot_ball.clone.res_y / 2), :] > 0)
            contact_delta = pressure_elements + 2 - self.rot_ball.clone.res_y
            contact_width += np.sign(contact_delta) * contact_width / 25

        self.rot_ball.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                                   self.init_force, contact_width=contact_width)
        self.stat_ball.make_slave_to(self.rot_ball)
        return init_displ

    def calc_contact_pressure(self, ui=None, res_dir=None):
        """Calculate contact pressure distribution between balls"""
        print_it('calculating 1 pressure distribution')

        init_displ = self.get_grid_size(ui, res_dir)
        [self.influ_mat_db_1] = load_influ_mat(ui, res_dir, 1)
        print_it('solving first half space', PrintOpts.lvl1.value)
        self.rot_ball.press, self.influ_mat_db_1 = \
            solve_half_space(self.rot_ball.profile, self.stat_ball.profile,
                             self.rot_ball.x_axis, self.rot_ball.y_axis,
                             self.rot_ball.res_x, self.rot_ball.res_y,
                             self.rot_ball.delta_x, self.rot_ball.delta_y,
                             self.rot_ball.e, self.stat_ball.e,
                             self.rot_ball.ny, self.stat_ball.ny,
                             self.rot_ball.norm_forces[0], res_dir,
                             init_displ=init_displ,
                             influ_mat_db=self.influ_mat_db_1)

        cache_influ_mat(ui, [self.influ_mat_db_1], res_dir)
        self.rot_ball.max_press = np.amax(self.rot_ball.press, axis=1)
        self.save_pressure_to_matlab(res_dir)

    def save_pressure_to_matlab(self, res_dir):
        dat_dict = dict(contact_pressure=self.rot_ball.press,
                        x_axis=self.rot_ball.x_axis,
                        y_axis=self.rot_ball.y_axis)
        save_to_matlab(dat_dict, res_dir, 'pressure-field')

    def calc_kinematics(self, rot_vel1, rot_vel2, ui=None, res_dir=None):
        """Calculate bearing kinematics based on rotational velocities of
        rotating balls"""
        print_it("calculating kinematics")
        self.effective_omega = rot_vel1 / 60
        self.raceway_vel = self.effective_omega * self.sliding_diam
        # neglecting sliding speed differences along ball contact area
        self.rel_vel = np.ones(self.rot_ball.res_x) * self.raceway_vel
        self.footpr_vel = self.raceway_vel
        self.overroll_t_incr = np.divide(self.rot_ball.delta_y, self.footpr_vel)
        self.press_zone_len = (self.rot_ball.press > 0).sum(
            1) * self.rot_ball.delta_y
        self.overroll_t = np.divide(self.press_zone_len, self.footpr_vel)
        self.no_overroll_t = np.divide((2 * math.pi * (
        self.sliding_diam / 2) - 3 * self.press_zone_len) / 3,
                                       self.footpr_vel)

    def calc_pv(self, ui=None, res_dir=None):
        """Calculate product of local maximum pressure and local maximum
        relative velocity"""
        print_it("calculating pv_rel")
        self.pv = np.multiply(self.rot_ball.max_press, self.rel_vel) / 1000

    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation in W per m^2"""
        print_it("calculating e_a,kin")
        pv_local = np.multiply(self.rot_ball.press.sum(1), self.rel_vel)
        self.e_akin = np.absolute(
            np.divide(np.multiply(pv_local, self.overroll_t_incr),
                      self.no_overroll_t)) / 1000

    def plot_it(self, ui=None, res_dir=None):
        """Orchestrate output plot generation"""
        print_it("plotting results and finishing up")
        plt_profile(self.rot_ball, PltOpts.DD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.rot_ball, PltOpts.DDD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.stat_ball, PltOpts.DD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.stat_ball, PltOpts.DDD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile_approx(res_dir, SubDir.profiles.value)
        plt_contact(self.rot_ball, self.stat_ball, PltOpts.DD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.rot_ball, self.stat_ball, PltOpts.DDD.value, res_dir,
                    SubDir.contacts.value)
        plt_3d(self.rot_ball.x_axis, self.rot_ball.y_axis, self.rot_ball.press,
               self.rot_ball.x_label,
               self.rot_ball.y_label, 'pressure in MPa', 'contact_pressure',
               res_dir, SubDir.pressures.value,
               'contact-pressure')
        plt_2d_scatt_line(self.rot_ball.x_axis, self.pv, self.rot_ball.x_axis,
                          self.pv, self.rot_ball.x_label,
                          'pv_rel in W $\mathregular{m^{-2}}$', 'pv_rel',
                          res_dir, SubDir.energy.value,
                          'pv_rel')
        plt_2d_scatt_line(self.stat_ball.x_axis, self.e_akin,
                          self.stat_ball.x_axis, self.e_akin,
                          self.stat_ball.x_label, 'a_akin in W $m^{-2}$',
                          'e_akin', res_dir,
                          SubDir.energy.value, 'e_akin')
        plt_2d_2y_ax(self.stat_ball.x_axis, self.e_akin, self.stat_ball.x_axis,
                     self.pv, self.stat_ball.x_label,
                     'e_akin vs pv_rel', 'e_akin in W $\mathregular{m^{-2}}$',
                     'pv_rel in W $\mathregular{m^{-2}}$',
                     res_dir, SubDir.energy.value, 'e_akin_vs_pv_rel')

    def generate_latex_output(self, calc_spec_tex_file_handle, sim, ui=None,
                              res_dir=None):
        """Generate calculation-specific part of the LaTeX output file"""
        average_pressure = np.mean(self.rot_ball.press[self.rot_ball.press > 0])
        numeric_output_data = [
            ('pressure, max.', to_preci(np.amax(self.rot_ball.press), 4),
             UnitTex.pressure.value, 'unverified'),
            ('pressure, av.', to_preci(average_pressure, 4),
             UnitTex.pressure.value, 'unverified'),
            ('e_a,kin, max.', to_preci(np.amax(self.e_akin), 4),
             UnitTex.eakin.value, 'unverified'),
            ('pv_rel, max.', to_preci(np.amax(self.pv), 4),
             UnitTex.pvrel.value, 'unverified'),
            ('contact area',
             to_preci(self.rot_ball.get_area(self.rot_ball.press), 4),
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
            TexTempl.FourBall.value, sim)
        with open(calc_spec_tex_file_handle, 'w') as f:
            f.write(template_calc_specific.render(latex_variables))

    def generate_latex_figures(self, ui=None, res_dir=None):
        """Generate calculation-specific figures for LaTeX report"""
        plt_contact(self.rot_ball, self.stat_ball, PltOpts.DDD.value, res_dir,
                    SubDir.tex_figs.value, 'contact1')

        plt_profile_approx(res_dir, SubDir.tex_figs.value)

        plt_3d(self.rot_ball.x_axis, self.rot_ball.y_axis, self.rot_ball.press,
               self.rot_ball.x_label,
               self.rot_ball.y_label, 'pressure in MPa', 'contact_pressure',
               res_dir, SubDir.tex_figs.value,
               'pressure1')

        plt_2d_2y_ax(self.stat_ball.x_axis, self.e_akin, self.stat_ball.x_axis,
                     self.pv, self.stat_ball.x_label,
                     'e_akin vs pv_rel', 'e_akin in W $\mathregular{m^{-2}}$',
                     'pv_rel in W $\mathregular{m^{-2}}$',
                     res_dir, SubDir.tex_figs.value, 'energy1')
