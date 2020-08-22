"""
Feature Model 2D
Main program
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from FeatMod2d_ops import width, height, res_x, res_z, num_ptcl
from FeatMod2d_mesh import MESHGRID
from FeatMod2d_ptcl import PARTICLE
from FeatMod2d_move import find_intersect_node
from FeatMod2d_plot import plot_mesh, plot_itsct

def hit_check(posn, mesh):
    int_x, int_z = np.divide(Arp.posn, 
                             np.array([res_x, res_z])).astype(int)
    return mesh.mat[int_x, int_z], (int_x, int_z)

threshold = 0.1
def update_surf(idx, mesh):
    if mesh.mat[idx] == 2:
            rnd = np.random.uniform(0.0, 1.0)
            if rnd < threshold:
                mesh.mat[idx] = 0


# create mesh
mesh = MESHGRID(width, height, res_x, res_z)
print(mesh)
mesh.mat_input()
mesh.find_surface()
mesh.plot()

record = [[] for i in range(num_ptcl)]

delta_L = min(res_x, res_z)
Arp = PARTICLE('Ar+', 'Ion',  32.0,     1)
for k in range(num_ptcl):
    Arp.init_posn(width, height)
    Arp.init_vels('Normal')
#    Arp.init_plot()
    record[k].append(Arp.posn.copy())
    
    for i in range(1000):
        Arp.move_ptcl(delta_L)
        Arp.posn[0] = Arp.posn[0] % mesh.width
        record[k].append(Arp.posn.copy())
        hit_mat, hit_idx = hit_check(Arp.posn, mesh)
        if hit_mat:
            break
    update_surf(hit_idx, mesh)
    record[k] = np.array(record[k]).T

#print(record[-1])

fig, ax = plt.subplots(1,1, figsize=(2,8),
                           constrained_layout=True)
colMap = cm.Accent
colMap.set_under(color='white')
x = np.linspace(0.0, mesh.width, mesh.nx) 
z = np.linspace(0.0, mesh.height, mesh.nz) 
X, Z = np.meshgrid(x, z) 
ax.contourf(X, Z, mesh.mat.T, cmap = colMap, vmin = 0.2, extend='both')
for i in range(20):
    ax.plot(record[num_ptcl-i-1][0,:], record[num_ptcl-i-1][1,:], 
            marker = 'o', markersize=1, linestyle='None' )
plt.show(fig)