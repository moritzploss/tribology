# Option for variable 'simulation_type':
# 1: cylindrical roller bearing
# 2:
# 3: cylindrical roller thrust bearing
# 4: ball on disk (currently not fully supported)
# 5: pin on disk
# 6: 4 ball
# 7: ball on three plates
# 8: ring on ring

# global simulation information
simulation_type = 6  # one of the above type numbers
simulation_name = 'FourBall'
auto_print = True  # True or False
auto_plot = True
auto_report = False  # reporting currently not supported, don't change

# global test setup / bearing information
tribo_system_name = 'foo'

# Rotating Ball (CB1)
e_cb1 = 210000  # young's modulus in MPa
ny_cb1 = 0.3  # poisson number in MPa
diameter_cb1 = 12.7  # in mm

# Stationary Balls (CB2)
e_cb2 = 210000  # young's modulus in MPa
ny_cb2 = 0.3  # poisson number in MPa
diameter_cb2 = 12.7  # in mm

# Loads
global_force = 50  # in N
rot_velocity1 = 50  # in rpm

# Mesh
res_x = 17  # data points along roller length
res_y = 15  # data points along roller width
