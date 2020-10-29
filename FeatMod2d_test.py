"""Feature Model 2D. Main program."""

# import numpy as np
from FeatMod2d_ops import (width, height, res_x, res_z, num_ptcl, ibc, 
                          threshold, max_rflct, idstrb, step_fac, max_step)
from FeatMod2d_mesh import MESHGRID
# from FeatMod2d_ptcl import PARTICLE
from Species import Arp
from FeatMod2d_mat import Si2d

# create mesh
mesh = MESHGRID(width, height, res_x, res_z)
print(mesh)
mesh.mat_input(Si2d)
mesh.find_surf()

delta_L = min(res_x, res_z)*step_fac

for k in range(num_ptcl):
    # print in process
    if (k + 1) % int(num_ptcl/5) == 0:
        print('%d particles are launched!' % (k+1))
        mesh.plot(dpi=300, fname='nptcl=%d.png' % (k+1))

    Arp.dead = 0
    Arp.init_posn(width, height)
    Arp.init_uvec(idstrb)
    # Arp.init_uvec(['Mono', 15.0])
    num_rflct = 0

    for i in range(max_step):
        # advance the ptcl by delta_L
        Arp.move_ptcl(delta_L)
        # periodic b.c. at left and right bdry
        Arp.bdry_check(mesh.width, mesh.height, ibc)
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
                mesh.update_surf(hit_idx)
                
        # check if the ptcl is dead
        if Arp.dead:
            break
