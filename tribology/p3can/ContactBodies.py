import math
import os

import numpy as np

from Constants import PrintOpts, Profiles
from artificial_surf import artificial_surf
from check_environment import print_import_error
from generate_profile import create_2d_profile2, revolve_profile2
from hertz_equations import get_diam_approx
from system_functions import print_it


class ContactBody:
    """Each tribosystems is composed of several contact bodies.
    The ContactBody class is the super class for all types of contact bodies,
    currently rings, balls, disks and flats"""

    # each contact body may use (and reuse) a separate
    # matlab engine for rough surface generation
    matlab_engine = None

    def __init__(self, name, e, ny):
        """All contact bodies have the following attributes"""
        self.name = name
        self.e = e
        self.ny = ny

        self.circle_radius = None
        self.clone = None
        self.contact_width = None
        self.diameter = None
        self.length = None
        self.normal_force = None
        self.path_profile = None
        self.profile = None
        self.rms_roughness = None
        self.rotational_speed = None
        self.roughness_mat = None
        self.true_rms_roughness = None
        self.type_profile = None

        self.res_x = None
        self.delta_x = None
        self.x_axis = None
        self.x_label = None
        self.x_profile = None
        self.file_x_axis = None
        self.file_x_profile = None
        self.r_hertz_x = None

        self.res_y = None
        self.delta_y = None
        self.y_axis = None
        self.y_label = None
        self.y_profile = None
        self.r_hertz_y = None

        self.res_z = None
        self.delta_z = None
        self.z_axis = None
        self.z_label = None
        self.z_profile = None

    def make_slave(self, contact_body):
        """Changes the contact body (slave) grid size and resolution to match
        that of another contact body (master)"""
        self.length = contact_body.length
        self.contact_width = contact_body.contact_width

        self.res_x = contact_body.res_x
        self.delta_x = contact_body.delta_x
        self.x_axis = contact_body.x_axis
        self.x_label = contact_body.x_label

        self.res_y = contact_body.res_y
        self.delta_y = contact_body.delta_y
        self.y_axis = contact_body.y_axis
        self.y_label = contact_body.y_label

        self.res_z = contact_body.res_z
        self.delta_z = contact_body.delta_z
        self.z_axis = contact_body.z_axis
        self.z_label = contact_body.z_label

    def get_hertz_radii(self):
        """Get an approximation of the hertz contact radii from the true
        contact body geometry"""
        self.r_hertz_x = get_diam_approx(self.x_axis, self.x_profile) / 2
        self.r_hertz_y = get_diam_approx(self.y_axis, self.y_profile) / 2

    def get_profile(self):
        """create a 2D profile, then revolve the 2D profile to obtain 3D
        profile"""
        self.x_axis, self.x_profile, self.file_x_axis, self.file_x_profile, self.delta_x = \
            create_2d_profile2(self.diameter, self.length, self.type_profile,
                               self.x_axis, self.res_x,
                               self.path_profile, self.circle_radius)
        self.profile, self.y_profile = \
            revolve_profile2(self.diameter, self.x_profile, self.y_axis,
                             self.res_x, self.res_y)

    def get_area(self, pressure_field):
        """Calculate true contact area based on number of non-zero pressure
        elements"""
        return len(
            pressure_field[pressure_field > 0]) * self.delta_x * self.delta_y

    def add_rms_roughness(self):
        """Add roughness to the contact body geometry. Currently, the roughness
        distribution is generated in matlab, hence, a matlab engine is required.
        A python implementation is on the TODO list"""
        if self.matlab_engine is None:
            try:
                print_it(
                    "trying to start matlab engine for {} "
                    "rough surface generation; this might take a while".
                    format(self.name), PrintOpts.lvl1.value)
                import matlab.engine
                self.matlab_engine = matlab.engine.start_matlab()
                self.matlab_engine.cd(os.getcwd(), nargout=0)
            except (ImportError, OSError):
                print_import_error('matlab engine for python',
                                   PrintOpts.lvl1.value)
                print_it(
                    "skipping roughness generation due to "
                    "missing matlab engine",
                    PrintOpts.lvl1.value)
                self.matlab_engine = 'no engine found'
        elif self.matlab_engine != 'no engine found':
            sigma = self.rms_roughness / 1000 / 1000
            h = 0.8
            m = self.res_y + 1
            lx = (self.contact_width + self.contact_width * self.res_y / (
            self.res_y + 1)) / 1000
            n = round(self.length / self.contact_width * m)
            qr = (2 * math.pi) / (self.contact_width / 20)
            z, pixel_width, qx, qy, cq, q, c, n, m = \
                artificial_surf(sigma, h, lx, m, n, qr, self.matlab_engine)

            extract_n = np.round(np.linspace(0, n - 1, self.res_x)).astype(int)
            self.roughness_mat = z[extract_n, 0:-1]
            self.true_rms_roughness = np.sqrt(
                np.mean(np.square(self.roughness_mat)))
            self.profile = np.add(self.profile, self.roughness_mat)
        else:
            pass


class Flat(ContactBody):
    """Contact body with completely flat surface and no specific physical size.
    This contact body has no 'make_profile' method since it is only supposed to
    be used as a slave to a Ball or Ring type contact body"""

    def __init__(self, name, e, ny, rms_rough=None):
        super().__init__(name, e, ny)
        self.type_profile = Profiles.NoProfile.value
        self.rms_roughness = rms_rough

    def make_slave_to(self, contact_body):
        """function call to 'make_slave' method of super class, plus some
        Flat-specific adjustments"""
        self.make_slave(contact_body)
        self.x_profile = np.zeros(self.res_x)
        self.y_profile = np.zeros(self.res_y)
        self.profile = np.zeros((self.res_x, self.res_y))
        self.get_hertz_radii()
        if self.rms_roughness is not None:
            self.add_rms_roughness()

    def simple_clone(self):
        """make a simple clone of contact body that is used for grid size
        determination"""
        self.clone = Flat('simple clone', self.e, self.ny)


class Ball(ContactBody):
    """Contact body with ball geometry"""

    def __init__(self, name, e, ny, diameter, rms_rough=None):
        super().__init__(name, e, ny)
        self.diameter = diameter
        self.type_profile = Profiles.Circle.value
        self.circle_radius = self.diameter / 2
        self.length = diameter
        self.r_hertz_y = diameter
        self.rms_roughness = rms_rough

    def make_profile(self, res_x, res_y, normal_force, contact_width=None):
        """mainly function call to 'get_profile' method of super class, plus
        some Ball-specific adjustments"""
        self.res_x = res_x
        self.res_y = res_y
        self.x_label = 'length in mm'
        self.y_label = 'width in mm'
        self.z_label = 'profile in mm'
        self.contact_width = \
            contact_width or get_contact_width(self.diameter, normal_force,
                                               self.ny, self.e, self.length)
        self.x_axis, self.delta_x = np.linspace(-self.contact_width / 2,
                                                self.contact_width / 2, res_x,
                                                retstep=True)
        self.y_axis, self.delta_y = np.linspace(-self.contact_width / 2,
                                                self.contact_width / 2, res_y,
                                                retstep=True)

        self.get_profile()
        self.get_hertz_radii()
        if self.rms_roughness is not None:
            self.add_rms_roughness()

    def make_slave_to(self, contact_body):
        """function call to 'make_slave' method of super class, plus some
        Ball-specific adjustments"""
        self.make_slave(contact_body)
        self.get_profile()
        self.get_hertz_radii()
        if self.rms_roughness is not None:
            self.add_rms_roughness()

    def simple_clone(self):
        """make a simple clone of contact body that is used for grid size
        determination"""
        self.clone = Ball('simple clone', self.e, self.ny, self.diameter)


class Disk(ContactBody):
    """Contact body with disk geometry. Difference to Flat-type contact body is
    that the physical size is explicitly
    defined here"""

    def __init__(self, name, e, ny, diameter=0, length=None, rot_vel=0,
                 rms_rough=None):
        super().__init__(name, e, ny)
        self.diameter = diameter
        self.rot_vel = rot_vel
        self.rms_roughness = rms_rough
        self.length = length

    def make_profile(self, res_x, res_y, normal_force, contact_width=None):
        """mainly function call to 'get_profile' method of super class, plus
        some Disk-specific adjustments"""
        self.res_x = res_x
        self.res_y = res_y
        self.x_label = 'length in mm'
        self.y_label = 'width in mm'
        self.z_label = 'profile in mm'
        self.contact_width = \
            contact_width or get_contact_width(self.diameter, normal_force,
                                               self.ny, self.e, self.length)
        self.x_axis, self.delta_x = np.linspace(-self.contact_width / 2,
                                                self.contact_width / 2, res_x,
                                                retstep=True)
        self.y_axis, self.delta_y = np.linspace(-self.contact_width / 2,
                                                self.contact_width / 2, res_y,
                                                retstep=True)

        self.x_profile = np.zeros(self.res_x)
        self.y_profile = np.zeros(self.res_y)
        self.profile = np.zeros((self.res_x, self.res_y))

        self.get_hertz_radii()
        if self.rms_roughness is not None:
            self.add_rms_roughness()

    def make_slave_to(self, contact_body):
        """function call to 'make_slave' method of super class, plus some
        Disk-specific adjustments"""
        self.make_slave(contact_body)
        self.x_profile = np.zeros(self.res_x)
        self.y_profile = np.zeros(self.res_y)
        self.profile = np.zeros((self.res_x, self.res_y))
        self.get_hertz_radii()
        if self.rms_roughness is not None:
            self.add_rms_roughness()

    def simple_clone(self):
        """make a simple clone of contact body that is used for grid size
        determination"""
        self.clone = Disk('simple clone', self.e, self.ny)


class Ring(ContactBody):
    """Ring type contact body. Can be used to generate ring, roller and
    arbitrarily-shaped contact bodies"""

    def __init__(self, name, e, ny, diameter, length, type_profile,
                 rms_rough=None, circle_radius=None,
                 path_profile=None):
        super().__init__(name, e, ny)
        self.type_profile = type_profile
        self.diameter = diameter
        self.length = length
        self.r_hertz_y = diameter
        self.psi = 0
        self.profile = None
        self.rms_roughness = rms_rough
        self.circle_radius = circle_radius
        self.path_profile = path_profile

    def make_profile(self, res_x, res_y, normal_force, contact_width=None,
                     contact_length=None):
        """mainly function call to 'get_profile' method of super class, plus
        some Ring-specific adjustments"""
        self.res_x = res_x
        self.res_y = res_y
        self.x_label = 'length in mm'
        self.y_label = 'width in mm'
        self.z_label = 'profile in mm'
        self.contact_width = \
            contact_width or get_contact_width(self.diameter, normal_force,
                                               self.ny, self.e, self.length)
        if contact_length is None:
            contact_length = self.length
        self.x_axis, self.delta_x = np.linspace(-contact_length / 2,
                                                contact_length / 2, res_x,
                                                retstep=True)
        self.y_axis, self.delta_y = np.linspace(-self.contact_width / 2,
                                                self.contact_width / 2, res_y,
                                                retstep=True)

        self.get_profile()
        self.get_hertz_radii()
        if self.rms_roughness is not None:
            self.add_rms_roughness()

    def make_slave_to(self, contact_body):
        """function call to 'make_slave' method of super class, plus some
        Ring-specific adjustments"""
        self.make_slave(contact_body)
        self.get_profile()
        self.get_hertz_radii()
        if self.rms_roughness is not None:
            self.add_rms_roughness()

    def simple_clone(self):
        """make a simple clone of contact body that is used for grid size
        determination"""
        self.clone = Ring('simple clone', self.e, self.ny, self.diameter,
                          self.length, Profiles.ISO.value)


def get_contact_width(diameter, normal_force, ny, e, length):
    """Very rough hertz contact width calculation for initial guess of grid
    size"""
    r_hertz = diameter / 2
    b = math.sqrt(
        1.5 * normal_force * r_hertz * (1 - math.pow(ny, 2)) / (2 * e * length))
    return 2 * b
