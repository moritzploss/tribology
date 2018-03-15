import math
import re

import numpy as np

from BluPrintTriboSys import TriboSys
from Constants import PltOpts, SubDir, TexTempl, UnitTex, PrintOpts, Unit, \
    PreSol
from cartesian_plot_functions import plt_contact, plt_profile, plt_3d, plt_2d, \
    plt_2d_scatt_line, plt_2d_2y_ax, \
    plt_profile_approx
from generate_latex_output import get_calc_specific_latex_template
from hertz_equations import hertz_displ
from influ_matrix_management import cache_influ_mat, load_influ_mat
from solve_half_space import solve_half_space, pre_solve_half_space
from system_functions import print_it, to_preci, save_to_matlab


class BallOnDisk(TriboSys):
    """Global tribosystem data"""

    def __init__(self, num_balls, sliding_diam, ball, disk, global_force,
                 setup_name='ball-on-disk'):
        super().__init__(setup_name, global_force, None)

        self.num_balls = num_balls
        self.sliding_diam = sliding_diam  # in mm
        self.ball = ball
        self.disk = disk
        self.rot_vel1 = None
        self.effective_rot_vel = None
        self.effective_omega = None
        self.raceway_vel = None
        self.ball.norm_forces = None
        self.ball.press = None
        self.ball.max_press = None
        self.ball.rel_vel = None
        self.pv = None
        self.init_force = None
        self.influ_mat_db_1 = None
        self.slip = None
        self.rel_vel = None
        if self.ball.diameter == self.ball.circle_radius:
            self.ball.is_ball = True
        else:
            self.ball.is_ball = False

    def calc_load_distribution(self, ui=None, res_dir=None):
        """Calculate load distribution"""
        print_it("calculating load distribution")
        self.ball.norm_forces = np.multiply(np.ones(self.num_balls),
                                            self.global_force / self.num_balls)
        self.init_force = self.global_force / self.num_balls

    def get_grid_size(self, ui=None, res_dir=None):
        """Determine grid size by running (quick) simulation with simplified
        contact bodies"""
        print_it('determining grid size', PrintOpts.lvl1.value)
        self.ball.simple_clone()
        self.ball.clone.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                                     self.init_force)
        self.disk.simple_clone()
        self.disk.clone.make_slave_to(self.ball.clone)

        init_displ = hertz_displ(self.ball.clone.e, self.disk.clone.e,
                                 self.ball.clone.ny, self.disk.clone.ny,
                                 self.ball.clone.r_hertz_x,
                                 self.ball.clone.r_hertz_y,
                                 self.disk.clone.r_hertz_x,
                                 self.disk.clone.r_hertz_y,
                                 self.ball.norm_forces[0])
        contact_delta = 1
        contact_width = 0.025
        while contact_delta != 0:
            self.ball.clone.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                                         self.init_force,
                                         contact_width=contact_width)
            self.disk.clone.make_slave_to(self.ball)

            pressure, init_displ = \
                pre_solve_half_space(self.ball.clone.profile,
                                     self.disk.clone.profile,
                                     self.ball.clone.x_axis,
                                     self.ball.clone.y_axis,
                                     self.ball.clone.res_x,
                                     self.ball.clone.res_y,
                                     self.ball.clone.delta_x,
                                     self.ball.clone.delta_y, self.ball.clone.e,
                                     self.disk.clone.e, self.ball.clone.ny,
                                     self.disk.clone.ny,
                                     self.ball.clone.norm_forces[0],
                                     init_displ=init_displ)

            pressure_elements = sum(
                pressure[math.floor(self.ball.clone.res_y / 2), :] > 0)
            contact_delta = pressure_elements + 2 - self.ball.clone.res_y
            contact_width += np.sign(contact_delta) * contact_width / 25

        self.ball.make_profile(PreSol.res_x.value, PreSol.res_y.value,
                               self.init_force, contact_width=contact_width)
        self.disk.make_slave_to(self.ball)
        return init_displ

    def calc_contact_pressure(self, ui=None, res_dir=None):
        """Calculate contact pressure distribution between ball and disk"""
        print_it('calculating 1 pressure distribution')

        init_displ = self.get_grid_size(ui, res_dir)
        [self.influ_mat_db_1] = load_influ_mat(ui, res_dir, 1)
        print_it('solving first half space', PrintOpts.lvl1.value)
        self.ball.press, self.influ_mat_db_1 = \
            solve_half_space(self.ball.profile, self.disk.profile,
                             self.ball.x_axis, self.ball.y_axis,
                             self.ball.res_x,
                             self.ball.res_y, self.ball.delta_x,
                             self.ball.delta_y, self.ball.e,
                             self.disk.e,
                             self.ball.ny, self.disk.ny,
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
        washer disks"""
        print_it("calculating kinematics")
        self.ball.rot_vel = rot_vel1
        self.disk.rot_vel = rot_vel2
        self.ball.omega = self.ball.rot_vel / 60
        self.disk.omega = self.disk.rot_vel / 60
        self.ball.vel = (
                        self.ball.diameter + self.ball.x_axis) * \
                        math.pi * self.ball.omega
        self.disk.vel = (
                        self.disk.diameter + self.ball.x_axis) * \
                        math.pi * self.disk.omega
        self.slip = (self.ball.vel - self.disk.vel) / self.ball.vel
        self.rel_vel = np.ones(self.ball.res_x) * (
        self.ball.vel - self.disk.vel)

        self.ball.footpr_vel = 2 * math.pi * self.ball.diameter / 2 * \
                               self.ball.omega
        self.disk.footpr_vel = 2 * math.pi * self.disk.diameter / 2 * \
                               self.disk.omega
        self.ball.overroll_t_incr = self.ball.delta_y / self.ball.footpr_vel
        self.disk.overroll_t_incr = self.disk.delta_y / self.disk.footpr_vel

        self.ball.overroll_t = np.divide(self.ball.press_zone_len,
                                         self.ball.footpr_vel)
        self.disk.overroll_t = np.divide(self.ball.press_zone_len,
                                         self.disk.footpr_vel)
        self.ball.no_overroll_t = np.divide(
            (2 * math.pi * (self.ball.diameter / 2) - self.num_balls *
             self.ball.press_zone_len) / self.num_balls, self.ball.footpr_vel)
        self.disk.no_overroll_t = np.divide(
            (2 * math.pi * (self.disk.diameter / 2) - self.ball.press_zone_len),
            self.disk.footpr_vel)

    def calc_pv(self, ui=None, res_dir=None):
        """Calculate product of local maximum pressure and local maximum
        relative velocity"""
        print_it("calculating pv_rel")
        self.pv = np.multiply(self.ball.max_press, self.rel_vel) / 1000

    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation in W per m^2"""
        print_it("calculating e_a,kin")
        pv_local = np.multiply(self.ball.press.sum(1), self.rel_vel)
        self.ball.e_akin = np.absolute(
            np.divide(np.multiply(pv_local, self.ball.overroll_t_incr),
                      self.ball.no_overroll_t)) / 1000
        self.disk.e_akin = np.absolute(
            np.divide(np.multiply(pv_local, self.disk.overroll_t_incr),
                      self.disk.no_overroll_t)) / 1000

    def plot_it(self, ui=None, res_dir=None):
        """Orchestrate output plot generation"""
        print_it("plotting results and finishing up")
        plt_profile(self.ball, PltOpts.DD.value, res_dir, SubDir.profiles.value)
        plt_profile(self.ball, PltOpts.DDD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.disk, PltOpts.DD.value, res_dir, SubDir.profiles.value)
        plt_profile(self.disk, PltOpts.DDD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile_approx(res_dir, SubDir.profiles.value)
        plt_contact(self.ball, self.disk, PltOpts.DD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.ball, self.disk, PltOpts.DDD.value, res_dir,
                    SubDir.contacts.value)

        plt_3d(self.ball.x_axis, self.ball.y_axis, self.ball.press,
               self.ball.x_label, self.ball.y_label,
               'pressure in MPa', 'contact_pressure_ball', res_dir,
               SubDir.pressures.value,
               'contact_pressure')
        plt_2d(self.ball.x_axis, self.rel_vel, self.ball.x_label, 'rel vel',
               'rel vel', res_dir,
               SubDir.energy.value, 'rel_vel')
        plt_2d(self.ball.x_axis, self.pv, self.ball.x_label,
               'pv_rel in {}'.format(Unit.pvrel.value), 'pv_rel',
               res_dir, SubDir.energy.value, 'pv_rel')
        plt_2d_scatt_line(self.ball.x_axis, self.ball.e_akin, self.ball.x_axis,
                          self.ball.e_akin,
                          self.ball.x_label,
                          'e_akin in {}'.format(Unit.eakin.value),
                          'e_akin ball',
                          res_dir, SubDir.energy.value, 'e-akin-ball')
        plt_2d_2y_ax(self.ball.x_axis, self.ball.e_akin, self.ball.x_axis,
                     self.pv, self.ball.x_label,
                     'e_akin vs pv_rel',
                     'e_akin in {}'.format(Unit.eakin.value),
                     'pv_rel in {}'.format(Unit.pvrel.value), res_dir,
                     SubDir.energy.value,
                     'e_akin_vs_pv_rel ball')

    def generate_latex_output(self, calc_spec_tex_file_handle, sim, ui=None,
                              res_dir=None):
        """Generate the calculation-specific part of the latex output file"""
        average_pressure = np.mean(self.ball.press[self.ball.press > 0])
        numeric_output_data = [
            ('pressure, max.', to_preci(np.amax(self.ball.press), 4),
             UnitTex.pressure.value, 'unverified'),
            ('pressure, av.', to_preci(average_pressure, 4),
             UnitTex.pressure.value, 'unverified'),
            ('e_a,kin ball, max.', to_preci(np.amax(self.ball.e_akin), 4),
             UnitTex.eakin.value, 'unverified'),
            ('e_a,kin disk, max.', to_preci(np.amax(self.disk.e_akin), 4),
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
                           'contact_plot1': '{}/contact1.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value),
                           'pressure_plot1': '{}/pressure1.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value),
                           'energy_plot1': '{}/energy1.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value),
                           'energy_plot2': '{}/energy2.png'.format(
                               SubDir.tex_figs_rel_to_tex_file.value)}
        template_calc_specific = get_calc_specific_latex_template(
            TexTempl.BallOnDisk.value, sim)
        with open(calc_spec_tex_file_handle, 'w') as f:
            f.write(template_calc_specific.render(latex_variables))

    def generate_latex_figures(self, ui=None, res_dir=None):
        """Generate calculation-specific part of the LaTeX output file"""
        plt_profile_approx(res_dir, SubDir.tex_figs.value)
        plt_contact(self.ball, self.disk, PltOpts.DDD.value, res_dir,
                    SubDir.tex_figs.value, 'contact1')
        plt_3d(self.ball.x_axis, self.ball.y_axis, self.ball.press,
               self.ball.x_label, self.ball.y_label,
               'pressure in {}'.format(UnitTex.pressure.value),
               'contact_pressure_ball', res_dir,
               SubDir.tex_figs.value, 'pressure1')
        plt_2d_2y_ax(self.ball.x_axis, self.ball.e_akin, self.ball.x_axis,
                     self.pv, self.ball.x_label,
                     'e_akin vs pv_rel',
                     'e_akin in {}'.format(UnitTex.eakin.value),
                     'pv_rel in {}'.format(UnitTex.pvrel.value), res_dir,
                     SubDir.tex_figs.value, 'energy1')
        plt_2d_2y_ax(self.disk.x_axis, self.disk.e_akin, self.disk.x_axis,
                     self.pv, self.disk.x_label,
                     'e_akin vs pv_rel',
                     'e_akin in {}'.format(UnitTex.eakin.value),
                     'pv_rel in {}'.format(UnitTex.pvrel.value), res_dir,
                     SubDir.tex_figs.value, 'energy2')
