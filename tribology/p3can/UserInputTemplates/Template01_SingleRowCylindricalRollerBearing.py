# Option for variable 'simulation_type':
# 1: cylindrical roller bearing
# 2:
# 3: cylindrical roller thrust bearing
# 4: ball on disk (currently not fully supported)
# 5: pin on disk
# 6: 4 ball
# 7: ball on three plates
# 8: ring on ring

# global simulation setup
simulation_type = 1  # one of the above type numbers
simulation_name = 'RadialBearing'
auto_print = True  # True or False
auto_plot = True
auto_report = False  # reporting currently not supported, don't change

# global test setup / bearing information
tribo_system_name = 'foo'
radial_clearance = 0.0184  # in mm
number_rollers = 13

#  Roller (CB1)
e_cb1 = 210000  # young's modulus in MPa
ny_cb1 = 0.3  # poisson number in MPa
diameter_cb1 = 9  # in mm
length_cb1 = 10  # in mm
type_profile_cb1 = 'ISO'  # 'None', 'ISO', 'Circle', 'File'
path_profile_cb1 = 'tribology/p3can/BearingProfiles/NU206-RE-1.txt'  # path to profile.txt file required if TypeProfile == 'File'
profile_radius_cb1 = 6.35  # input required if TypeProfile == 'Circle'
roller_slip = "tribology/p3can/RollerSlip/testslip.txt"

#  Inner Ring (CB2)
e_cb2 = 210000  # young's modulus in MPa
ny_cb2 = 0.3  # poisson number in MPa
diameter_cb2 = 37.5  # in mm
type_profile_cb2 = 'None'  # 'None', 'File'
profile_radius_cb2 = 20
length_cb2 = 10
path_profile_cb2 = "tribology/p3can/BearingProfiles/NU206-IR-2.txt"  # path to profile.txt file required if TypeProfile == 'File'

# Outer Ring (CB3)
e_cb3 = 210000  # young's modulus in MPa
ny_cb3 = 0.3  # poisson number in MPa
diameter_cb3 = 55.5  # in mm
type_profile_cb3 = 'None'  # 'None', 'File'
length_cb3 = 10
profile_radius_cb3 = 5
path_profile_cb3 = "tribology/p3can/BearingProfiles/NU2208-IR-1.txt"  # path to profile.txt file required if TypeProfile == 'File'

# Loads
global_force = 3000  # in N
rot_velocity1 = 0  # in rpm
rot_velocity2 = 300  # in rpm

# Mesh
res_x = 17  # data points along roller length
res_y = 15  # data points along roller width
res_pol = 13  # data points along bearing circumference
