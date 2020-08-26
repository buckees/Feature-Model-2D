"""
Feature Model 2D
Main program
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from FeatMod2d_ops import width, height, res_x, res_z, num_ptcl, threshold
from FeatMod2d_mesh import MESHGRID
from FeatMod2d_ptcl import PARTICLE
from FeatMod2d_rflct import REFLECT

def rotate_random(uvec):
    theta = np.random.uniform(-math.pi/4.0, +math.pi/4.0)
    theta = theta + math.pi
    x1, z1 = uvec
    x2 = math.cos(theta)*x1 - math.sin(theta)*z1
    z2 = math.sin(theta)*x1 + math.cos(theta)*z1
    return np.array([x2, z2])

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

    imove_ptcl, ireaction, iremove = 1, 0, 0
    num_reflect = 0
    while imove_ptcl == 1 and num_reflect < 5:
        for i in range(1000):
            Arp.move_ptcl(delta_L)
            Arp.posn[0] = Arp.posn[0] % mesh.width
            record[k].append(Arp.posn.copy())
            if Arp.posn[1] >= mesh.height:
                iremove = 1
                imove_ptcl = 0
                break
            
            hit_mat, hit_idx = mesh.hit_check(Arp.posn)
            if hit_mat:
                break
            
        if iremove == 1:
            break
        
        if not hit_mat:
            num_reflect += 1
            continue
        # at this position, th ptcl hits a mat
        # decide wehter a reflection or reaction
        
        rnd = np.random.uniform(0.0, 1.0)
        mat_name = mesh.mater[hit_mat]
        rflct = REFLECT(Arp.name, mat_name, 1.0)
        prob = rflct.calc_prob()
        if rnd < prob:
    #        call reflection
            u1 = Arp.uvec
            Arp.uvec = rotate_random(Arp.uvec)
            num_reflect += 1
            u2 = Arp.uvec
#            angle = np.arccos(np.clip(np.dot(-u1, u2), -1, 1))
#            angle = angle/math.pi*180.0
#            print(angle)
        else:
            imove_ptcl = 0
            ireaction = 1
            
    if ireaction == 1:
        mesh.update_surf(hit_idx, threshold)
    
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
for i in range(10):
    ax.plot(record[num_ptcl-i-1][0,:], record[num_ptcl-i-1][1,:], 
            marker = 'o', markersize=1, linestyle='None' )
plt.show(fig)