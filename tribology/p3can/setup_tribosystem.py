import math

from BluPrintAxialThrustBear import AxialThrustBear
from BluPrintBallOnThreePlates import BallOnThreePlates
from BluPrintCylRollBear import CylRollBear
from BluPrintFourBall import FourBall
from BluPrintPinOnDisk import PinOnDisk
from BluPrintRingOnRing import RingOnRing
from Constants import Profiles
from ContactBodies import Ring, Flat, Disk, Ball
from system_functions import exit_program, in_d


def setup_cyl_rol_bearing(ui):
    """Generate roller contact body (master) as well as inner and outer ring
    (slaves)"""
    roller = Ring('roller', ui['e_cb1'], ui['ny_cb1'], ui['diameter_cb1'],
                  ui['length_cb1'], ui['type_profile_cb1'],
                  rms_rough=in_d('rms_roughness_cb1', ui),
                  circle_radius=in_d('profile_radius_cb1', ui),
                  path_profile=in_d('path_profile_cb1', ui))
    roller.make_profile(ui['res_x'], ui['res_y'], ui['global_force'], 1)
    ring1 = Ring('inner-ring', ui['e_cb2'], ui['ny_cb2'], ui['diameter_cb2'],
                 ui['length_cb2'], ui['type_profile_cb2'],
                 rms_rough=in_d('rms_roughness_cb2', ui),
                 circle_radius=in_d('profile_radius_cb2', ui),
                 path_profile=in_d('path_profile_cb2', ui))
    ring1.make_slave_to(roller)
    ring2 = Ring('outer-ring', ui['e_cb3'], ui['ny_cb3'], -ui['diameter_cb3'],
                 ui['length_cb3'], ui['type_profile_cb3'],
                 rms_rough=in_d('rms_roughness_cb2', ui),
                 circle_radius=in_d('profile_radius_cb3', ui),
                 path_profile=in_d('path_profile_cb3', ui))
    ring2.make_slave_to(roller)
    tribo_system = CylRollBear(ui['number_rollers'], roller, ring1, ring2,
                               ui['global_force'], ui['radial_clearance'],
                               in_d('res_pol', ui), in_d('roller_slip', ui),
                               in_d('tribo_system_name', ui))
    return roller, ring1, ring2, tribo_system


def setup_deep_gro_ball_bearing(ui):
    """Not currently implemented"""
    exit_program("simulation_type '2' currently undefined")


def setup_cyl_rol_thrust_bearing(ui):
    """Generate roller contact body (master) as well as two flats (slave)"""
    roller = Ring('roller', ui['e_cb1'], ui['ny_cb1'], ui['diameter_cb1'],
                  ui['length_cb1'], ui['type_profile_cb1'],
                  rms_rough=in_d('rms_roughness_cb1', ui),
                  circle_radius=in_d('profile_radius_cb1', ui),
                  path_profile=in_d('path_profile_cb1', ui))
    roller.make_profile(ui['res_x'], ui['res_y'],
                        ui['global_force'] / ui['number_rollers'], 1)
    ring1 = Flat('shaft-washer', ui['e_cb2'], ui['ny_cb2'],
                 rms_rough=in_d('rms_roughness_cb2', ui))
    ring1.make_slave_to(roller)
    ring2 = Flat('housing-washer', ui['e_cb3'], ui['ny_cb3'],
                 rms_rough=in_d('rms_roughness_cb3', ui))
    ring2.make_slave_to(roller)
    tribo_system = AxialThrustBear(ui['number_rollers'], ui['mean_diameter'],
                                   roller, ring1, ring2,
                                   ui['global_force'],
                                   in_d('tribo_system_name', ui))
    return roller, ring1, ring2, tribo_system


def setup_ball_on_disk(ui):
    """Not currently implemented"""
    exit_program("simulation_type '4' currently not implemented completely.")
    """
    ball = Ball('pin', ui['e_cb1'], ui['ny_cb1'], ui['diameter_cb1'],
                rms_rough=in_d('rms_roughness_cb1', ui))
    ball.make_profile(ui['res_x'], ui['res_y'],
                      ui['global_force'] / ui['number_balls'])
    disk = Disk('disk', ui['e_cb2'], ui['ny_cb2'], diameter=ui['diameter_cb2'],
                rot_vel='rot_velocity2',
                rms_rough=in_d('rms_roughness_cb2', ui))
    disk.make_slave_to(ball)
    tribo_system = BallOnDisk(ui['number_balls'], ui['sliding_diameter'], ball,
                              disk, ui['global_force'],
                              in_d('tribo_system_name', ui))   
    return ball, disk, tribo_system
    """
    pass


def setup_pin_on_disk(ui):
    """Generate pin contact body (master; either ball or ring,
    depending on user input), then generate disk (slave)"""
    if ui['type_profile_cb1'] == Profiles.Circle.value:
        if ui['diameter_cb1'] == 2 * ui['profile_radius_cb1']:
            pin = Ball('pin', ui['e_cb1'], ui['ny_cb1'], ui['diameter_cb1'],
                       rms_rough=in_d('rms_roughness_cb1', ui))
            pin.make_profile(ui['res_x'], ui['res_y'],
                             ui['global_force'] / ui['number_pins'])
        else:
            pin = Ring('pin', ui['e_cb1'], ui['ny_cb1'], ui['diameter_cb1'],
                       ui['length_cb1'], ui['type_profile_cb1'],
                       rms_rough=in_d('rms_roughness_cb1', ui),
                       circle_radius=in_d('profile_radius_cb1', ui),
                       path_profile=in_d('path_profile_cb1', ui))
            pin.make_profile(ui['res_x'], ui['res_y'],
                             ui['global_force'] / ui['number_pins'], 1)
    else:
        pin = Ring('pin', ui['e_cb1'], ui['ny_cb1'], ui['diameter_cb1'],
                   ui['length_cb1'], ui['type_profile_cb1'],
                   rms_rough=in_d('rms_roughness_cb1', ui),
                   circle_radius=in_d('profile_radius_cb1', ui),
                   path_profile=in_d('path_profile_cb1', ui))
        pin.make_profile(ui['res_x'], ui['res_y'],
                         ui['global_force'] / ui['number_pins'], 1)
    disk = Disk('disk', ui['e_cb2'], ui['ny_cb2'],
                rms_rough=in_d('rms_roughness_cb2', ui))
    disk.make_slave_to(pin)
    tribo_system = PinOnDisk(ui['number_pins'], ui['sliding_diameter'], pin,
                             disk, ui['global_force'],
                             in_d('tribo_system_name', ui))
    return pin, disk, tribo_system


def setup_four_ball(ui):
    """Generate rotating ball (master) and stationary ball (slave) contact
    bodies"""
    rot_ball = Ball('rotating-ball', ui['e_cb1'], ui['ny_cb1'],
                    ui['diameter_cb1'],
                    rms_rough=in_d('rms_roughness_cb1', ui))
    rot_ball.make_profile(ui['res_x'], ui['res_y'],
                          ui['global_force'] / math.sqrt(6))
    stat_ball = Ball('stationary-ball', ui['e_cb2'], ui['ny_cb2'],
                     ui['diameter_cb2'],
                     rms_rough=in_d('rms_roughness_cb2', ui))
    stat_ball.make_slave_to(rot_ball)
    tribo_system = FourBall(rot_ball, stat_ball, ui['global_force'],
                            in_d('tribo_system_name', ui))
    return rot_ball, stat_ball, tribo_system


def setup_ball_on_three_plates(ui):
    """Generate ball (master) and flat (slave) contact bodies"""
    ball = Ball('ball', ui['e_cb1'], ui['ny_cb1'], ui['diameter_cb1'],
                rms_rough=in_d('rms_roughness_cb1', ui))
    ball.make_profile(ui['res_x'], ui['res_y'],
                      ui['global_force'] * math.sqrt(2) / 3)
    plate = Flat('plate', ui['e_cb2'], ui['ny_cb2'],
                 rms_rough=in_d('rms_roughness_cb2', ui))
    plate.make_slave_to(ball)
    tribo_system = BallOnThreePlates(ball, plate, ui['global_force'],
                                     in_d('plate_angle', ui, math.pi / 2),
                                     in_d('tribo_system_name', ui))
    return ball, plate, tribo_system


def setup_ring_on_ring(ui):
    """"Generate central ring contact body (sun) as well as planetary ring(s),
    then make planets slave to sun"""
    sun = Ring('sun', ui['e_cb1'], ui['ny_cb1'], ui['diameter_cb1'],
               ui['length_cb1'], ui['type_profile_cb1'],
               rms_rough=in_d('rms_roughness_cb1', ui),
               circle_radius=in_d('profile_radius_cb1', ui),
               path_profile=in_d('path_profile_cb1', ui))
    sun.make_profile(ui['res_x'], ui['res_y'],
                     ui['global_force'] / ui['number_planets'], 1)
    planet = Ring('planet', ui['e_cb2'], ui['ny_cb2'], ui['diameter_cb2'],
                  ui['length_cb2'], ui['type_profile_cb2'],
                  rms_rough=in_d('rms_roughness_cb2', ui),
                  circle_radius=in_d('profile_radius_cb2', ui),
                  path_profile=in_d('path_profile_cb2', ui))
    planet.make_slave_to(sun)
    tribo_system = RingOnRing(ui['number_planets'], sun, planet,
                              ui['global_force'], in_d('tribo_system_name', ui))
    return sun, planet, tribo_system
