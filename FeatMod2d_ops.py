"""Feature Model 2D. Operation parameters."""

# mesh informaiton
width, height = 100.0, 260.0  # nm
res_x, res_z = 1.0, 1.0  # nm

# particle information
num_ptcl = 100000

# reflection information
max_rflct = -1

# reaction informaiton
threshold = 1.0

# init uvec distribution
idstrb = ['Uniform3D', -45.0, 45.0]

# boundary condition
ibc = 'periodic'

# step length factor
step_fac = 0.5

# max steps for a single particle
max_step = 1000

# surf norm calc
surf_norm_range = 4
surf_norm_mode = 'Sum Vector'