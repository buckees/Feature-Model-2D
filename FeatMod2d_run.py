# -*- coding: utf-8 -*-
"""
Feature Model 2D
run file
"""
import matplotlib.pyplot as plt
import numpy as np
import math

from FeatMod2d_ops import *
from FeatMod2d_geo import *

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
        mesh[ncoord[0]:ncoord[1], ncoord[2]:ncoord[3]] = mater[1]
        
fig, axes = plt.subplots(1,2, figsize=(4,12),
                           constrained_layout=True)
    
ax1 = axes[0]
ax1.contourf(mesh.T)

mesh_surf, surf = find_surface(mesh)

ax2 = axes[1]
ax2.contourf(mesh_surf.T)

#fig2, axes2 = plt.subplots(1,2, figsize=(4,4),
#                           constrained_layout=True)
#ax1 = axes2[0]
#ax1.plot(surf)

