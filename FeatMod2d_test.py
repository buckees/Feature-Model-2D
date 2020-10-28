"""Feature Model 2D. Main program."""

import numpy as np
from FeatMod2d_ops import width, height, res_x, res_z, num_ptcl, \
                          threshold, max_rflct
from FeatMod2d_mesh import MESHGRID
# from FeatMod2d_ptcl import PARTICLE
from Species import Arp
from FeatMod2d_mat import Si2d

# create mesh
res_x, res_z = 0.5, 0.5
mesh = MESHGRID(width, height, res_x, res_z)
print(mesh)
mesh.mat_input(Si2d)

delta_L = min(res_x, res_z)*0.5

for k in range(num_ptcl):
    # print in process
    if (k + 1) % int(num_ptcl/5) == 0:
        print('%d particles are launched!' % (k+1))
        mesh.plot(dpi=300, fname='nptcl=%d.png' % (k+1))

    Arp.dead = 0
    Arp.init_posn(width, height)
    Arp.init_uvec(['Uniform', -15.0, 15.0])
    num_rflct = 0

    for i in range(10000):
        # advance the ptcl by delta_L
        Arp.move_ptcl(delta_L)
        # periodic b.c. at left and right bdry
        Arp.bdry_check(mesh.width, mesh.height, 'periodic')
        # check if the ptcl is dead
        if Arp.dead:
            break
        hit_mat, hit_idx = mesh.hit_check(Arp.posn)
        if hit_mat:
            mat_name = mesh.mater[hit_mat]
            Arp.dead = 1
            # calc surf norm
            if mat_name == 'Si':
                # now ireact = 1
                mesh.mat[hit_idx] = 0
                
        # check if the ptcl is dead
        if Arp.dead:
            break
