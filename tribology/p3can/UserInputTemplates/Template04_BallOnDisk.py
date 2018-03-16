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
simulation_type = 4  # one of the above types
simulation_name = 'BallOnDisk'
auto_print = True  # True or False
auto_plot = True
auto_report = False  # reporting currently not supported, don't change

# global test setup / bearing information
tribo_system_name = 'foo'
number_balls = 1
sliding_diameter = 50  # in mm, sliding diameter on disc

# Pin (CB1)
e_cb1 = 210000  # young's modulus in MPa
ny_cb1 = 0.3  # poisson number in MPa
diameter_cb1 = 12.7

# Disk (CB2)
e_cb2 = 210000  # young's modulus in MPa
ny_cb2 = 0.3  # poisson number in MPa
diameter_cb2 = 50

# Loads
global_force = 50  # in N
rot_velocity1 = 50  # in rpm
rot_velocity2 = 10

# Mesh
res_x = 33  # data points along roller length
res_y = 31  # data points along roller width
