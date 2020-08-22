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

# create mesh
mesh = MESHGRID(width, height, res_x, res_z)
print(mesh)
mesh.mat_input()
mesh.find_surface()
mesh.plot()

Arp = PARTICLE('Ar+', 'Ion',  32.0,     1, num_ptcl)
print(Arp)
Arp.init_posn(width, height)
Arp.init_vels()
Arp.init_plot()


def hit_check(posn, mesh):
    int_x, int_z = np.floor(Arp.posn)[:,0].astype(int)
    ibdry = 0
    if not (0 < int_x < mesh.nx-1):
        ibdy = 1
    return mesh.mat[int_x, int_z], (int_x, int_z), ibdry

delta_L = min(res_x, res_z)
for i in range(100):
    Arp.move_ptcl(delta_L)
    hit_mat, hit_idx, ibdry = hit_check(Arp.posn, mesh)
    if hit_mat or ibdry: break
Arp.init_plot()


