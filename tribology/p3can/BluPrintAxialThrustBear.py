import copy
import math
import re

import numpy as np

from BluPrintRollBear import AxialRollBear
from Constants import PltOpts, SubDir, MasterDir, UnitTex, Unit, PrintOpts, \
    PreSol
from cartesian_plot_functions import plt_profile, plt_contact, plt_3d, \
    plt_2d_scatt_line, plt_2d_2y_ax, plt_profile_approx
from hertz_equations import hertz_displ
from influ_matrix_management import cache_influ_mat, load_influ_mat
from solve_half_space import solve_half_space, pre_solve_half_space
from system_functions import print_it, to_preci


class AxialThrustBear(AxialRollBear):
    """Global bearing data"""

    def __init__(self, num_rollers, mean_diameter, roller, ring1, ring2,
                 init_force, bearing_name='bearing'):
        super().__init__(bearing_name, num_rollers, roller, ring1, ring2,
                         init_force, mean_diameter)

        self.press_zone_len = None
        self.overrolling_duration = None
        self.overroll_t_incr = None
        self.non_overroll_t = None
        self.e_akin = None
        self.uniq_norm_forces = None
        self.eff_rot_velocity = None
        self.influ_mat_db_1 = None

    def calc_load_distribution(self, ui=None, res_dir=None):
        """Calculate load distribution"""
        print_it("calculating 1 load distribution")
        self.roller.norm_forces = \
            np.multiply(np.ones(self.num_rollers),
                        self.global_force / self.num_rollers)
        self.uniq_norm_forces = self.roller.norm_forces[0]

    def get_grid_size(self, ui=None, res_dir=None):
        """Determine grid size by running (quick) simulation with simplified
        contact bodies"""
        print_it('determining grid size', PrintOpts.lvl1.value)
        self.roller.simple_clone()
        self.roller.clone.make_profile(PreSol.res_x.value,
                                       PreSol.res_y.value, self.init_force)
        self.ring1.simple_clone()
        self.ring1.clone.make_slave_to(self.roller.clone)

        init_displ = \
            hertz_displ(self.roller.clone.e, self.ring1.clone.e,
                        self.roller.clone.ny, self.ring1.clone.ny,
                        self.roller.clone.r_hertz_x,
                        self.roller.clone.r_hertz_y,
                        self.ring1.clone.r_hertz_x,
                        self.ring1.clone.r_hertz_y, self.roller.norm_forces[0])
        contact_delta = 1
        contact_width = 0.05
        while contact_delta != 0:
            self.roller.clone.make_profile(self.roller.clone.res_x,
                                           self.roller.clone.res_y,
                                           self.init_force,
                                           contact_width=contact_width)
            self.ring1.clone.make_slave_to(self.roller.clone)

            pressure, init_displ = \
                pre_solve_half_space(self.roller.clone.profile,
                                     self.ring1.clone.profile,
                                     self.roller.clone.x_axis,
                                     self.roller.clone.y_axis,
                                     self.roller.clone.res_x,
                                     self.roller.clone.res_y,
                                     self.roller.clone.delta_x,
                                     self.roller.clone.delta_y,
                                     self.roller.clone.e,
                                     self.ring1.clone.e,
                                     self.roller.clone.ny,
                                     self.ring1.clone.ny,
                                     self.roller.norm_forces[0],
                                     init_displ=init_displ)

            pressure_elements = sum(
                pressure[math.floor(self.roller.clone.res_y / 2), :] > 0)
            contact_delta = pressure_elements + 2 - self.roller.clone.res_y
            contact_width += np.sign(contact_delta) * contact_width / 25

        self.roller.make_profile(self.roller.res_x, self.roller.res_y,
                                 self.init_force, contact_width=contact_width)
        self.ring1.make_slave_to(self.roller)
        self.ring2.make_slave_to(self.roller)
        return init_displ

    def calc_contact_pressure(self, ui=None, res_dir=None):
        """Calculate contact pressure distribution between washer(s) and
        roller"""
        print_it('calculating 1 pressure distribution')
        init_displ = self.get_grid_size(ui, res_dir)
        [self.influ_mat_db_1] = load_influ_mat(ui, res_dir, 1)
        print_it('solving first half space', PrintOpts.lvl1.value)
        self.roller.press, self.influ_mat_db_1 = \
            solve_half_space(self.roller.profile, self.ring1.profile,
                             self.roller.x_axis, self.roller.y_axis,
                             self.roller.res_x, self.roller.res_y,
                             self.roller.delta_x, self.roller.delta_y,
                             self.roller.e, self.ring1.e, self.roller.ny,
                             self.ring1.ny, self.roller.norm_forces[0],
                             res_dir, init_displ=init_displ,
                             influ_mat_db=self.influ_mat_db_1)
        self.roller.max_press = np.amax(self.roller.press, axis=1)
        self.ring1.press = copy.copy(self.roller.press)
        self.ring1.max_press = copy.copy(self.roller.max_press)
        cache_influ_mat(ui, [self.influ_mat_db_1], res_dir)

    def calc_kinematics(self, rot_vel1, rot_vel2, ui=None, res_dir=None):
        """Calculate bearing kinematics based on rotational velocities of
        washer disks"""
        print_it("calculating kinematics")
        self.rot_vel1 = rot_vel1
        self.rot_vel2 = rot_vel2
        self.eff_rot_velocity = self.rot_vel1 - self.rot_vel2
        self.eff_omega = self.eff_rot_velocity / 60
        self.eff_omega_cage = self.eff_omega / 2
        self.roller.eff_omega = self.mean_diameter / self.roller.diameter * \
                                self.eff_omega_cage
        self.footpr_vel = 2 * math.pi * \
                          (self.mean_diameter / 2 + self.roller.x_axis) * \
                          self.eff_omega / 2

        self.overroll_t_incr = self.roller.delta_y / self.footpr_vel
        self.press_zone_len = (self.ring1.press > 0).sum(1) * self.ring1.delta_y
        self.overrolling_duration = np.divide(self.press_zone_len,
                                              self.footpr_vel)
        self.non_overroll_t = np.divide(
            (2 * math.pi * (self.mean_diameter / 2 + self.roller.x_axis) -
             self.num_rollers * self.press_zone_len) / self.num_rollers,
            self.footpr_vel)

        self.ring1.eff_raceway_vel = 2 * math.pi * (
        self.mean_diameter / 2 + self.roller.x_axis) * self.eff_omega
        self.ring2.eff_raceway_vel = 0
        self.roller.actual_raceway_vel = 2 * math.pi * (
        self.mean_diameter / 2 + self.roller.x_axis) * \
                                 self.eff_omega_cage + math.pi * \
                                 self.roller.diameter * \
                                 self.roller.eff_omega
        self.roller.eff_raceway_vel = self.roller.actual_raceway_vel
        self.roller.rel_vel = math.pi * self.roller.diameter * \
                              self.roller.eff_omega - 2 * math.pi * \
                              (self.mean_diameter / 2 + self.roller.x_axis) * \
                              self.eff_omega_cage
        self.roller.slip = np.divide(self.roller.rel_vel, 2 * math.pi *
                                     (self.mean_diameter / 2 +
                                      self.roller.x_axis) *
                                     self.eff_omega_cage)

    def calc_pv(self, ui=None, res_dir=None):
        """Calculate product of local maximum pressure and local maximum
        relative velocity"""
        print_it("calculating pv_rel")
        self.pv = np.absolute(np.multiply(self.ring1.max_press,
                                          self.roller.rel_vel)) / 1000

    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation in W per m^2"""
        print_it("calculating e_a,kin")
        pv_local = np.multiply(self.ring1.press.sum(1), self.roller.rel_vel)
        self.e_akin = np.absolute(np.divide(np.multiply(pv_local,
                                                        self.overroll_t_incr),
                                            self.non_overroll_t)) / 1000

    def plot_it(self, ui=None, res_dir=None):
        """Orchestrate output plot generation"""
        print_it("plotting results and finishing up")
        self.plt_profiles_and_contacts(res_dir)
        self.plt_pol_slip(res_dir)
        self.plt_pres(res_dir)
        self.plt_pv(res_dir)
        self.plt_eakin(res_dir)

    def plt_profiles_and_contacts(self, res_dir):
        """Generate plots of profiles and contacts"""
        plt_profile(self.roller, PltOpts.DD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.roller, PltOpts.DDD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.ring1, PltOpts.DD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.ring1, PltOpts.DDD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.ring2, PltOpts.DD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile(self.ring2, PltOpts.DDD.value, res_dir,
                    SubDir.profiles.value)
        plt_profile_approx(res_dir, SubDir.profiles.value)
        plt_contact(self.roller, self.ring1, PltOpts.DD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.roller, self.ring1, PltOpts.DDD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.roller, self.ring2, PltOpts.DD.value, res_dir,
                    SubDir.contacts.value)
        plt_contact(self.roller, self.ring2, PltOpts.DDD.value, res_dir,
                    SubDir.contacts.value)

    def plt_pol_slip(self, res_dir):
        """Plot roller slip"""
        plt_2d_scatt_line(self.roller.x_axis, self.roller.slip,
                          self.roller.x_axis, self.roller.slip,
                          self.roller.x_label, "roller slip in percent",
                          "roller slip", res_dir,
                          SubDir.energy.value, "roller slip")

    def plt_pres(self, res_dir):
        """Plot contact pressure distribution"""
        plt_3d(self.roller.x_axis, self.roller.y_axis, self.ring1.press,
               self.roller.x_label, self.roller.y_label, 'pressure in MPa',
               'contact pressure', res_dir, SubDir.pressures.value,
               'contact-pressure')

    def plt_pv(self, res_dir):
        """Plot pv"""
        plt_2d_scatt_line(self.roller.x_axis, self.pv, self.roller.x_axis,
                          self.pv, self.roller.x_label,
                          'pv_rel in {}'.format(Unit.pvrel.value),
                          'pv_rel', res_dir, SubDir.energy.value, 'pv_rel')

    def plt_eakin(self, res_dir):
        """Plot kinetic friction energy accumulation"""
        plt_2d_scatt_line(self.roller.x_axis, self.e_akin, self.roller.x_axis,
                          self.e_akin, self.roller.x_label,
                          'e_akin in W $m^{-2}$', 'e_akin', res_dir,
                          SubDir.energy.value, 'e_akin')
        plt_2d_2y_ax(self.roller.x_axis, self.e_akin, self.roller.x_axis,
                     self.pv, self.roller.x_label, 'e_akin vs pv_rel',
                     'e_akin in {}'.format(Unit.eakin.value),
                     'pv_rel in {}'.format(Unit.pvrel.value), res_dir,
                     SubDir.energy.value, 'e_akin_vs_pv_rel')

    def generate_latex_output(self, calc_spec_tex_file_handle, sim, ui=None,
                              res_dir=None):
        """Generate calculation-specific part of the LaTeX output file"""
        mean_press = np.mean(self.ring1.press[:, :][self.ring1.press[:, :] > 0])
        numeric_output_data = [
            ('pv_rel, max.', to_preci(np.amax(self.pv), 4),
             UnitTex.pvrel.value, 'unverified'),
            ('e_akin, max.', to_preci(np.amax(self.e_akin), 4),
             UnitTex.eakin.value, 'unverified'),
            ('pressure, max.', to_preci(np.amax(self.ring1.press[:, :]), 4),
             UnitTex.pressure.value, 'unverified'),
            ('pressure, av.', to_preci(mean_press, 4),
             UnitTex.pressure.value, 'unverified'),
            ('slip, max.', to_preci(np.amax(self.roller.slip) * 100, 4),
             '\%', 'unverified'),
            ('normal force, max.', to_preci(max(self.roller.norm_forces), 4),
             UnitTex.force.value, 'unverified'),
            ('contact area', to_preci(self.roller.get_area(self.ring1.press),
                                      4),
             UnitTex.area.value, 'unverified')]

        scale_factor_ir = 0.5 + 0.1 * (self.rot_vel1 != 0)
        scale_factor_or = 0.5 + 0.1 * (self.rot_vel2 != 0)

        table_calc_summary = []
        for key, value, unit, status in sorted(numeric_output_data):
            table_calc_summary.append(
                (re.sub('_', '\_', key), value, unit, status)
            )

        latex_variables = {'contact_plot1': 'figures{}contact1.png'.format('/'),
                           'pressure_plot1':
                               'figures{}pressure1.png'.format('/'),
                           'slip_plot1': 'figures{}slip1.png'.format('/'),
                           'energy_plot1': 'figures{}energy1.png'.format('/'),
                           'scale_factor_ir': scale_factor_ir,
                           'scale_factor_or': scale_factor_or,
                           'table_calc_summary': table_calc_summary}

        template_calc_specific = sim.latex_jinja_env.get_template(
            MasterDir.latex_templates.value +
            '/' +
            'template-axial-thrust-bearing.tex')
        with open(calc_spec_tex_file_handle, 'w') as f:
            f.write(template_calc_specific.render(latex_variables))

    def generate_latex_figures(self, ui=None, res_dir=None):
        """Generate calculation-specific figures for LaTeX report"""
        plt_profile_approx(res_dir, SubDir.tex_figs.value)
        plt_contact(self.roller, self.ring1, PltOpts.DDD.value, res_dir,
                    SubDir.tex_figs.value, 'contact1')
        plt_3d(self.roller.x_axis, self.roller.y_axis, self.ring1.press[:, :],
               self.roller.x_label, self.roller.y_label,
               'pressure in MPa', 'pressure at highest normal force', res_dir,
               SubDir.tex_figs.value, 'pressure1')
        plt_2d_scatt_line(self.roller.x_axis, self.roller.slip,
                          self.roller.x_axis, self.roller.slip,
                          self.roller.x_label, "roller slip in percent (0-1)",
                          "roller slip", res_dir,
                          SubDir.tex_figs.value, "slip1")
        plt_2d_2y_ax(self.roller.x_axis, self.e_akin, self.roller.x_axis,
                     self.pv, self.roller.x_label, 'e_akin vs pv_rel',
                     'e_akin in {}'.format(Unit.eakin.value),
                     'pv_rel in {}'.format(Unit.pvrel.value), res_dir,
                     SubDir.tex_figs.value, 'energy1')
