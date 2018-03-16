# Option for variable 'simulation_type':
# 1: cylindrical roller bearing
# 2:
# 3: cylindrical roller thrust bearing
# 4: ball on disk
# 5: pin on disk
# 6: 4 ball
# 7: ball on three plates
# 8: ring on ring

# global simulation setup
simulation_type = 3  # one of the above type numbers
simulation_name = 'Sample-Simulation'
auto_print = True  # True or False
auto_plot = True
auto_report = False  # reporting currently not supported, don't change

# global test setup / bearing information
tribo_system_name = 'bla'
radial_clearance = 0.0184  # in mm
number_balls = 1
number_rollers = 13
number_pins = 1
number_planets = 1
sliding_diameter = 42
mean_diameter = 77.5

#  Contact Body 1 (CB1)
e_cb1 = 210000  # young's modulus in MPa
ny_cb1 = 0.3  # poisson number in MPa
diameter_cb1 = 9 # in mm
length_cb1 = 10  # in mm
type_profile_cb1 = 'ISO'  # 'None', 'ISO', 'Circle', 'File'
path_profile_cb1 = 'tribology/p3can/BearingProfiles/NU206-RE-1.txt'  # path to profile.txt file required if TypeProfile == 'File'
profile_radius_cb1 = 6.35  # input required if TypeProfile == 'Circle'
roller_slip = "tribology/p3can/RollerSlip/testslip.txt"

#  Contact Body 2 (CB2)
e_cb2 = 210000  # young's modulus in MPa
ny_cb2 = 0.3  # poisson number in MPa
diameter_cb2 = 37.5  # in mm
type_profile_cb2 = 'None'  # 'None', 'File'
profile_radius_cb2 = 20
length_cb2 = 10
path_profile_cb2 = "tribology/p3can/BearingProfiles/NU206-IR-2.txt"  # path to profile.txt file required if TypeProfile == 'File'

# Contact Body 3 (CB3)
e_cb3 = 210000  # young's modulus in MPa
ny_cb3 = 0.3  # poisson number in MPa
diameter_cb3 = 55.5  # in mm
type_profile_cb3 = 'None'  # 'None', 'File'
length_cb3 = 10
profile_radius_cb3 = 5
path_profile_cb3 = ''+''  # path to profile.txt file required if TypeProfile == 'File'

# Loads
global_force = 30000  # in N
rot_velocity1 = 300  # in rpm
rot_velocity2 = 0  # in rpm

# Mesh
res_x = 41  # data points along roller length
res_y = 41  # data points along roller width
res_pol = 15  # data points along bearing circumference
