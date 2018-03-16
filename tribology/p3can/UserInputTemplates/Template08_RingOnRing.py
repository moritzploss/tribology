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
simulation_type = 8  # one of the above type numbers
simulation_name = 'RingOnRing'
auto_print = True  # True or False
auto_plot = True
auto_report = False  # reporting currently not supported, don't change

# global test setup / bearing information
tribo_system_name = 'foo'
number_planets = 2

# Sun (CB1)
e_cb1 = 210000  # young's modulus in MPa
ny_cb1 = 0.3  # poisson number in MPa
diameter_cb1 = 37.5 # in mm
length_cb1 = 10  # in mm
type_profile_cb1 = 'ISO'  # 'None', 'ISO', 'Circle', 'File'
path_profile_cb1 = 'tribology/p3can/BearingProfiles/NU206-RE-1.txt'  # path to profile.txt file required if TypeProfile == 'File'
profile_radius_cb1 = 6.35  # input required if TypeProfile == 'Circle'

# Planet (CB2)
e_cb2 = 210000  # young's modulus in MPa
ny_cb2 = 0.3  # poisson number in MPa
diameter_cb2 = 9  # in mm
type_profile_cb2 = 'ISO'  # 'None', 'File'
profile_radius_cb2 = 20
length_cb2 = 10
path_profile_cb2 = "tribology/p3can/BearingProfiles/NU206-IR-2.txt"  # path to profile.txt file required if TypeProfile == 'File'

# Loads
global_force = 300  # in N
rot_velocity1 = 300  # in rpm
rot_velocity2 = 2  # in rpm

# Mesh
res_x = 33  # data points along roller length
res_y = 31  # data points along roller width
