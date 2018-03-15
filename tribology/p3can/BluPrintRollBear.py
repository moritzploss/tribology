from abc import ABCMeta, abstractmethod

from BluPrintTriboSys import TriboSys


class RollBear(TriboSys):
    # TODO: add class attribute describtion
    """

    Attributes:
        wheels: An integer representing the number of wheels the vehicle has.

    """

    __metaclass__ = ABCMeta

    def __init__(self, name, number_rollers, roller, inner_ring, outer_ring,
                 global_force, init_force):
        super().__init__(name, global_force, init_force)
        self.roller = roller
        self.ring1 = inner_ring
        self.ring2 = outer_ring
        self.num_rollers = number_rollers
        self.rot_vel1 = 0
        self.rot_vel2 = 0

        self.mean_diam = None
        self.eff_rot_vel = None
        self.eff_omega = None
        self.eff_omega_cage = None

        self.footpr_vel = None
        self.actual_raceway_vel_roller = None

    @abstractmethod
    def calc_load_distribution(self, ui=None, res_dir=None):
        """Calculate the load distribution within the tribosystem"""
        pass

    @abstractmethod
    def calc_contact_pressure(self, ui=None, res_dir=None):
        """Calculate the contact pressure for each contact in the tribosystem"""
        pass

    @abstractmethod
    def calc_kinematics(self, rot_velocity, rot_velocity2, ui=None,
                        res_dir=None):
        """Calculate the kinematics of the tribosystem"""
        pass

    @abstractmethod
    def calc_pv(self, ui=None, res_dir=None):
        """"Calculate the product of the local maximum pressure and the
        relative velocity between ring and roller"""
        pass

    @abstractmethod
    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation"""
        pass

    @abstractmethod
    def plot_it(self, ui=None, res_dir=None):
        """Plot relevant calculation results"""
        pass

    @abstractmethod
    def generate_latex_output(self, calc_spec_tex_file_handle, simulation,
                              ui=None, res_dir=None):
        pass

    @abstractmethod
    def generate_latex_figures(self, ui=None, res_dir=None):
        pass


class AxialRollBear(RollBear):
    """Global bearing data"""

    def __init__(self, name, number_rollers, roller, shaft_ring, housing_ring,
                 global_force, mean_diameter):
        super().__init__(name, number_rollers, roller, shaft_ring, housing_ring,
                         global_force,
                         global_force / number_rollers)

        self.mean_diameter = mean_diameter
        self.pv = None

    @abstractmethod
    def calc_load_distribution(self, ui=None, res_dir=None):
        """Calculate the load distribution within the tribosystem"""
        pass

    @abstractmethod
    def calc_contact_pressure(self, ui=None, res_dir=None):
        """Calculate the contact pressure for each contact in the tribosystem"""
        pass

    @abstractmethod
    def calc_kinematics(self, rot_velocity, rot_velocity2, ui=None,
                        res_dir=None):
        """Calculate the kinematics of the tribosystem"""
        pass

    @abstractmethod
    def calc_pv(self, ui=None, res_dir=None):
        """"Calculate the product of the local maximum pressure and the relative
         velocity between ring and roller"""
        pass

    @abstractmethod
    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation"""
        pass

    @abstractmethod
    def plot_it(self, ui=None, res_dir=None):
        """Plot relevant calculation results"""
        pass

    @abstractmethod
    def generate_latex_output(self, calc_spec_tex_file_handle, simulation,
                              ui=None, res_dir=None):
        pass

    @abstractmethod
    def generate_latex_figures(self, ui=None, res_dir=None):
        pass


class RadialRollBear(RollBear):
    """Global bearing data"""

    def __init__(self, name, number_rollers, roller, inner_ring, outer_ring,
                 global_force, radial_clearance, res_pol,
                 path_roller_slip):
        super().__init__(name, number_rollers, roller, inner_ring, outer_ring,
                         global_force, global_force)

        self.mean_diameter = (self.ring1.diameter - self.ring2.diameter) / 2
        self.rad_clear = radial_clearance
        self.pv_ring1 = None
        self.pv_ring2 = None
        self.path_roller_slip = path_roller_slip
        self.res_pol = res_pol
        self.pol_ax = None
        self.dd_phi = None
        self.d_phi = None
        self.num_roller_pos = None
        self.phi_mat = None
        self.phi = None

    @abstractmethod
    def calc_load_distribution(self, ui=None, res_dir=None):
        """Calculate the load distribution within the tribosystem"""
        pass

    @abstractmethod
    def calc_contact_pressure(self, ui=None, res_dir=None):
        """Calculate the contact pressure for each contact in the tribosystem"""
        pass

    @abstractmethod
    def calc_kinematics(self, rot_velocity, rot_velocity2, ui=None,
                        res_dir=None):
        """Calculate the kinematics of the tribosystem"""
        pass

    @abstractmethod
    def calc_pv(self, ui=None, res_dir=None):
        """"Calculate the product of the local maximum pressure and the relative
        velocity between ring and roller"""
        pass

    @abstractmethod
    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation"""
        pass

    @abstractmethod
    def calc_e_akin(self, ui=None, res_dir=None):
        """"Calculate the kinetic friction energy accumulation"""
        pass

    @abstractmethod
    def plot_it(self, ui=None, res_dir=None):
        """Plot relevant calculation results"""
        pass

    @abstractmethod
    def generate_latex_output(self, calc_spec_tex_file_handle, simulation,
                              ui=None, res_dir=None):
        """Generate LaTeX output"""
        pass

    @abstractmethod
    def generate_latex_figures(self, ui=None, res_dir=None):
        """Generate LaTeX figures"""
        pass
