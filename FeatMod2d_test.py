"""Feature Model 2D. Main program."""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from copy import copy, deepcopy
import numpy as np

# import numpy as np
from FeatMod2d_ops import (width, height, res_x, res_z, num_ptcl, ibc, 
                          threshold, max_rflct, idstrb, step_fac, max_step,
                          num_plot, surf_norm_range, surf_norm_mode)
from FeatMod2d_mesh import MESHGRID
# from FeatMod2d_ptcl import PARTICLE
from Species import Arp
from FeatMod2d_mat import Si2d_trench
from FeatMod2d_rflct import REFLECT

# create mesh
mesh = MESHGRID(width, height, res_x, res_z)
print(mesh)
mesh.mat_input(Si2d_trench)
mesh.find_surf()

delta_L = min(res_x, res_z)*step_fac
Arp_rflct = REFLECT()


# record for diagnostics
rec_traj = []

for k in range(num_ptcl):
    # print in process
    if (k + 1) % int(num_ptcl/num_plot) == 0:
        print('%d particles are launched!' % (k+1))
        mesh.plot(dpi=300, fname='nptcl=%d.png' % (k+1))
        mesh.plot_surf(surf_norm_range=surf_norm_range, 
                       surf_norm_mode=surf_norm_mode, 
                       dpi=300, fname='surf_nptcl=%d.png' % (k+1))
        mesh.find_float_cell(idiag=1)

    Arp.dead = 0
    Arp.init_posn(width, height)
    Arp.init_uvec(idstrb)
    num_rflct = 0
    
    # record the particle trajectory
    rec_traj.append([])
    rec_traj[-1].append(Arp.posn.copy())

    for i in range(max_step):
        # advance the ptcl by delta_L
        Arp.move_ptcl(delta_L)
        # periodic b.c. at left and right bdry
        Arp.bdry_check(mesh.width, mesh.height, ibc)
        # check if the ptcl is dead
        if Arp.dead:
            rec_traj[-1].append(Arp.posn.copy())
            break
        hit_mat, hit_idx = mesh.hit_check(Arp.posn)
        if hit_mat:
            rec_traj[-1].append(Arp.posn.copy())
            mat_name = mesh.mater[hit_mat]
            Arp_rflct.svec, Arp_rflct.stheta = \
                mesh.calc_surf_norm(hit_idx, radius=surf_norm_range, 
                                    imode=surf_norm_mode)
            # calc surf norm
            if mat_name == 'Si':
                # now ireact = 1
                mesh.mat[hit_idx] = 0
                mesh.update_surf(hit_idx)
                Arp.dead = 1
            else:
                if num_rflct > max_rflct:
                    Arp.dead = 1
                    rec_traj[-1].append(Arp.posn.copy())
                    break
                # Test diffusive rflct only
                Arp.uvec = Arp_rflct.diff_rflct()
                # Arp.uvec = Arp_rflct.spec_rflct(Arp.uvec)
                # move the ptcl by 10 steps until it gets out of the mat
                for ii in range(10):
                    Arp.move_ptcl(delta_L)
                    Arp.bdry_check(mesh.width, mesh.height, ibc)
                    hit_mat, hit_idx = mesh.hit_check(Arp.posn)
                    if not hit_mat:
                        break
                if hit_mat:
                    Arp.dead = 1
                    rec_traj[-1].append(Arp.posn.copy())
                    break
                num_rflct += 1
                
        # check if the ptcl is dead
        if Arp.dead:
            rec_traj[-1].append(Arp.posn.copy())
            break
        
    rec_traj[-1] = np.array(rec_traj[-1]).T

colMap = copy(cm.Accent)
colMap.set_under(color='white')

def plot_traj(ax, traj):
    for temp_ptcl in traj[-50::]:
        # print(temp_ptcl)
        ax.plot(temp_ptcl[0, :], temp_ptcl[1, :],
                marker='o', markersize=0.3, linestyle='-', linewidth=0.1)


fig, axes = plt.subplots(1, 2, figsize=(12, 8),
                          constrained_layout=True)

ax = axes[0]
ax.contourf(mesh.x, mesh.z, mesh.mat, cmap=colMap, vmin=0.2, extend='both')
ax.set_xlim(0.0, mesh.width)
ax.set_ylim(0.0, mesh.height)
plot_traj(ax, rec_traj)

ax = axes[1]
ax.scatter(mesh.x, mesh.z, c=mesh.mat, s=1, cmap=colMap, vmin=0.2)
ax.set_xlim(0.0, mesh.width)
ax.set_ylim(0.0, mesh.height)
plot_traj(ax, rec_traj)

fig.savefig('ptcl_traj.png', dpi=300)
