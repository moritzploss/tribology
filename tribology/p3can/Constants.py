import os
from enum import Enum
from matplotlib import cm


class SignCheck(Enum):
    """Options for signs of numeric variables"""
    no_check = 0
    pos = 1
    neg = -1


class OddEvenCheck(Enum):
    """Options for odd-even check of numeric variables"""
    Odd = 1
    Even = 0
    Undefined = None


class CMaps(Enum):
    """Colormap options"""
    default = cm.get_cmap('PuBu')
    markers = default
    contacts = default


class VarType(Enum):
    """Variable type options"""
    string = 0
    integer = 1
    real_number = 2
    bool = 3


class PrintOpts(Enum):
    """Print levels used to structure the console output"""
    lvl0 = ""
    lvl1 = "\t"
    lvl2 = "\t\t"


class SimType(Enum):
    """Types of tribosystems that can be simulated using p3can"""
    cyl_rol_bearing = 1
    deep_gro_ball_bearing = 2
    cyl_rol_thrust_bearing = 3
    ball_on_disk = 4
    pin_on_disk = 5
    four_ball = 6
    ball_on_three_plates = 7
    ring_on_ring = 8


class Profiles(Enum):
    """Profile types for contact bodies"""
    ISO = 'ISO'
    Circle = 'Circle'
    File = 'File'
    NoProfile = 'None'


class SubDir(Enum):
    """Subdirectories of the root directory that are used to save calculation
    outputs"""
    contacts = "plots-contact"
    profiles = "plots-profile"
    pressures = "plots-pressure"
    energy = "plots-energy-figures"
    load_distr = "plots-load-distribution"
    matlab = "matlab-files"
    latex_files = "latex-files"
    tex_figs_rel_to_tex_file = "figures"
    tex_figs = "{}{}{}".format(latex_files, os.sep, tex_figs_rel_to_tex_file)
    pub_figs = "publication-figures"
    np_db = "numpy-databases"
    infl_mat_db_folder = 'cache'


class MasterDir(Enum):
    """Names of top level root directories"""
    latex_templates = os.sep.join(os.path.realpath(__file__)
                                  .split(os.sep)[:-1]) \
                      + os.sep + "LatexTemplates"
    results = "Results"


class NpDBs(Enum):
    """Constants related to influence matrix storage and cache size"""
    infl_mat_db = 'infl-mat-db.pickle'  # name of influence matrix database
    max_cache_size = 1024  # maximum size of influence matrix cache in megabyte


class Unit(Enum):
    """Units of physical quantities for python/matplotlib output"""
    press = "MPa"
    force = "N"
    eakin = "W $\mathregular{m^{-2}}$"
    pvrel = "W $\mathregular{m^{-2}}$"
    area = "W $\mathregular{m^{-2}}$"


class UnitTex(Enum):
    """Units of physical quantities for LaTeX output"""
    pressure = "MPa"
    force = "N"
    eakin = "W m$^{-2}$"
    pvrel = "W m$^{-2}$"
    area = "mm$^{2}$"


class TexTempl(Enum):
    """Names of LaTeX templates"""
    CylindricalRollerBearing = "template-cylindrical-roller-bearing.tex"
    CylindricalRollerThrustBearing = "template-cylindrical-roller-" \
                                     "thrust-bearing.tex"
    PinOnDisk = "template-pin-on-disk.tex"
    FourBall = "template-four-ball.tex"
    BallOn3Plates = "template-ball-on-three-plates.tex"
    RingOnRing = "template-ring-on-ring.tex"
    BallOnDisk = "template-ball-on-disk.tex"
    Generic = 'p3can-calculation-report-template.tex'


class PreSol(Enum):
    """Grid resolution used for the pre-solver when determining a grid size for
     the actual half space calculation"""
    res_x = 13
    res_y = 15


class AutoReportOptions(Enum):
    """Options for automatic LaTeX report generation"""
    true = True
    false = False


class PltOpts(Enum):
    """Options for output plots"""
    DD = '2D'
    DDD = '3D'
    res = 96     # default output plot resolution in dpi
    msize = 75  # default marker size in scatter plots
    azim = 120  # default spacial orientation of non-standard 3D plots


class VarDicts(Enum):
    """Required and optional variables for ALL simulation types.
    Simulation-specific required and optional variables
    are handled elsewhere"""
    required_user_input = {'res_x': VarType.integer,
                           'res_y': VarType.integer,
                           'e_cb1': VarType.real_number,
                           'ny_cb1': VarType.real_number,
                           'diameter_cb1': VarType.real_number,
                           'e_cb2': VarType.real_number,
                           'ny_cb2': VarType.real_number,
                           'global_force': VarType.real_number,
                           'rot_velocity1': VarType.real_number,
                           'simulation_name': VarType.string,
                           'simulation_type': VarType.integer}
    optional_user_input = {'tribo_system_name': VarType.string,
                           'auto_print': VarType.bool,
                           'auto_plot': VarType.bool,
                           'auto_report': VarType.bool,
                           'rms_roughness_cb1': VarType.real_number,
                           'rms_roughness_cb2': VarType.real_number}


class VarProps(Enum):
    """Properties of variables that apply independent of simulation type"""
    positive_vars = ['res_x', 'res_y', 'res_z', 'e_cb1', 'ny_cb1',
                     'diameter_cb1', 'length_cb1', 'e_cb2', 'ny_cb2',
                     'diameter_cb2', 'length_cb2', 'e_cb3', 'ny_cb3',
                     'diameter_cb3', 'length_cb3', 'global_force', 'res_pol',
                     'number_rollers', 'number_pins', 'number_planets',
                     'sliding_diameter', 'mean_diameter', 'plate_angle',
                     'rms_roughness_cb1', 'rms_roughness_cb2',
                     'rms_roughness_cb3']
    negative_vars = []
    alphanumeric_vars = ['simulation_name', 'tribo_system_name']
    odd_vars = ['res_x', 'res_y', 'res_z']
    bool_vars = ['auto_print', 'auto_plot', 'auto_report']
    int_vars = ['number_pins', 'number_rollers', 'number_planets', 'res_x',
                'res_y', 'res_z', 'res_pol']
