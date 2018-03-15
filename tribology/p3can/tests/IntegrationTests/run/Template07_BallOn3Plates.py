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
simulation_type = 7  # one of the above type numbers
simulation_name = 'BallOnThreePlates'
auto_print = True  # True or False
auto_plot = True
auto_report = True

# global test setup / bearing information
tribo_system_name = 'bla'
plate_angle = 45  # angle in degree

# Rotating Ball (CB1)
e_cb1 = 210000  # young's modulus in MPa
ny_cb1 = 0.3  # poisson number in MPa
diameter_cb1 = 12.7  # in mm

# Stationary Plates (CB2)
e_cb2 = 210000  # young's modulus in MPa
ny_cb2 = 0.3  # poisson number in MPa

# Loads
global_force = 50  # in N
rot_velocity1 = 50  # in rpm

# Mesh
res_x = 17  # data points along roller length
res_y = 15  # data points along roller width
