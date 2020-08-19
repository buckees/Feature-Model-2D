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
from FeatMod2d_plot import plot_traj

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
        print(ncoord)
        mesh[ncoord[0]:ncoord[1], ncoord[2]:ncoord
             [3]] = mater[1]
        
fig, axes = plt.subplots(1,3, figsize=(6,12),
                           constrained_layout=True)
    
colMap = cm.Accent
colMap.set_under(color='white')

ax0 = axes[0]
ax0.contourf(mesh.T, cmap = colMap,vmin = 0.2, extend='both')

mesh_surf, node_surf = find_surface(mesh)

ax1 = axes[1]
ax1.contourf(mesh_surf.T)
ax2 = axes[2]
ax2.plot(node_surf[0, :], node_surf[1, :], 'o-')

plt.show(fig)

#fig2, axes2 = plt.subplots(1,2, figsize=(4,4),
#                           constrained_layout=True)
#ax1 = axes2[0]
#ax1.plot(surf)

# initiate particles
posx = np.random.uniform(0.0, 1.0, size=num_ptcl)*domain_x
posz = np.ones(num_ptcl)*domain_z
posn_init = np.array([posx, posz])
#theta = np.random.uniform(0.0, 1.0, size=num_ptcl)*180.0
ux, uz = random_2d_vector(num_ptcl)

fig2, axes2 = plt.subplots(1,2, figsize=(8,3),
                           constrained_layout=True)
ax0 = axes2[0]
ax0.plot(ux, uz, 'o')
ax1 = axes2[1]
ax1.quiver(0.0, 0.0, ux, uz)
plt.show(fig2)

# (x = posx + ux*L, z = posz + uz*L)
# line intersects with x=0
# write in form of ax + bz + c = 0
# uz*x - ux*z -(uz*posx - ux*posz) = 0

#line_a, line_b, line_c = uz, -ux, -(uz*posx-ux*posz)
line = np.array([uz, -ux, -(uz*posx-ux*posz)])
node_itsct = find_intersect_node(line, node_surf)

fig3, axes3 = plt.subplots(1,3, figsize=(6,12),
                           constrained_layout=True)

ax0 = axes3[0]
ax0.contourf(mesh.T, cmap = colMap,vmin = 0.2, extend='both')
ax0.plot(node_itsct[0, :], node_itsct[1, :], 'ro')

for i in range(num_ptcl):
    point1 = posn_init[:, i]
    point2 = node_itsct[:, i]
    plot_traj(ax0, point1, point2)

plt.show(fig3)
