# -*- coding: utf-8 -*-
"""
Feature Model 2D
run file
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from FeatMod2d_ops import *
from FeatMod2d_geo import *
from FeatMod2d_init import random_2d_vector
from FeatMod2d_move import find_intersect_node
from FeatMod2d_plot import plot_mesh, plot_itsct

# create geometry and mesh
nx, nz = (math.ceil(domain_x/resolution), math.ceil(domain_z/resolution))
print('nx = %d' % nx, ', nz = %d' % nz)

# create a 2d array as mesh
mesh = np.zeros((nx, nz))
# assign materials to the mesh
for mat in materials:
    mater = mat[0]
    print(mater)
    itype = mat[1]
    if itype == 'rect':
        coord = mat[2]
        ncoord = rect_conv(coord)
        mesh[ncoord[0]:ncoord[1], ncoord[2]:ncoord
             [3]] = mater[1]

mesh_surf, node_surf = find_surface(mesh)

# plot initial mesh and surface
plot_mesh(mesh, node_surf)

# initiate particles
posx = np.random.uniform(0.0, 1.0, size=num_ptcl)*domain_x
posz = np.ones(num_ptcl)*domain_z
posn_init = np.array([posx, posz])
#theta = np.random.uniform(0.0, 1.0, size=num_ptcl)*180.0
ux, uz = random_2d_vector(num_ptcl)
uvec = np.array([ux, uz])
# (x = posx + ux*L, z = posz + uz*L)
L = 100
node_traj = posn_init + uvec*L

# (x = posx + ux*L, z = posz + uz*L)
# line intersects with x=0
# write in form of ax + bz + c = 0
# uz*x - ux*z -(uz*posx - ux*posz) = 0

#line_a, line_b, line_c = uz, -ux, -(uz*posx-ux*posz)
line = np.array([uz, -ux, -(uz*posx-ux*posz)])
node_itsct = find_intersect_node(line, node_surf)

plot_itsct(mesh, posn_init, node_traj, node_itsct)