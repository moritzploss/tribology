from abc import ABCMeta, abstractmethod


class TriboSys(object):
    """All simulations rely on the definition of an object of type TriboSys,
    which contains all physical parameters describing the tribosystem,
    including contact bodies."""

    __metaclass__ = ABCMeta

    def __init__(self, name, global_force, init_force):
        self.name = name
        self.init_force = init_force
        self.global_force = global_force

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
        """Calculate the kinematics of the tribosystem"""
        pass

    @abstractmethod
    def calc_e_akin(self, ui=None, res_dir=None):
        """Calculate the kinematics of the tribosystem"""
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
