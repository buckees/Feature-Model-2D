"""Feature Model 2D. Main program."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from copy import copy, deepcopy

from FeatMod2d_ops import (width, height, res_x, res_z, num_ptcl, ibc, 
                          threshold, max_rflct, idstrb, step_fac, max_step,
                          num_plot, surf_norm_range, surf_norm_mode)
from FeatMod2d_mesh import MESHGRID
from FeatMod2d_ptcl import PARTICLE
from Species import Arp
from FeatMod2d_rflct import REFLECT
from FeatMod2d_mat import Si2d

# create mesh
mesh = MESHGRID(width, height, res_x, res_z)
print(mesh)
mesh.add_mat(Si2d)
# mesh.plot()

delta_L = min(res_x, res_z)*step_fac
# delta_L = min(res_x, res_z)

# species information is imported from species
# Initialize the PARTICLE() object
ptcl = PARTICLE(**Arp)
ptcl_rflct = REFLECT()

# init diagnostics
rec_traj, rec_surf, rec_mesh = [], [], []

for k in range(num_ptcl):
    if (k + 1) % int(num_ptcl/num_plot) == 0:
        print('%d particles are launched!' % (k+1))
        mesh.plot(dpi=300, fname='nptcl=%d.png' % (k+1))
        # rec_mesh.append(deepcopy(mesh.mat))
    ptcl.dead = 0
    ptcl.init_posn(width, height)
    # record initial position
    if k > num_ptcl - 20:
        rec_traj.append([])
        rec_traj[-1].append(ptcl.posn.copy())
    ptcl.init_uvec(idstrb)

    num_rflct = 0
#    while imove_ptcl == 1 and num_rflct < 5:
    for i in range(max_step):
        # advance the ptcl by delta_L
        ptcl.move_ptcl(delta_L)
        # periodic b.c. at left and right bdry
        ptcl.bdry_check(mesh.width, mesh.height, ibc)
        # check if the ptcl is dead
        if ptcl.dead:
            # record ptcl posn when dead
            if k > num_ptcl - 20:
                rec_traj[-1].append(ptcl.posn.copy())
            break
        hit_mat, hit_idx = mesh.hit_check(ptcl.posn)
        if hit_mat:
            # record the hit point
            if k > num_ptcl - 20:
                rec_traj[-1].append(ptcl.posn.copy())
            # at this position, th ptcl hits a mat
            mat_name = mesh.mater[hit_mat]
            # calc surf norm
            ptcl_rflct.svec, ptcl_rflct.stheta = \
                mesh.calc_surf_norm(hit_idx, radius=surf_norm_range, 
                                    imode=surf_norm_mode)
            # decide wehter a reflection or reaction
            rnd = np.random.uniform(0.0, 1.0)
            rflct = REFLECT(ptcl.name, mat_name, 1.0)
            prob = rflct.calc_prob()
            if rnd < prob:
                # check max rflct
                if num_rflct > max_rflct:
                    ptcl.dead = 1
                    break
                
                ptcl.uvec = ptcl_rflct.rflct(ptcl.uvec)
                num_rflct += 1
            else:
                # now ireact = 1
                mesh.update_mat(hit_idx, threshold)
                mesh.find_float_cell()
                ptcl.dead = 1
        # check if the ptcl is dead
        if ptcl.dead:
            # record ptcl posn when dead
            if k > num_ptcl - 20:
                rec_traj[-1].append(ptcl.posn.copy())
            break

    if k > num_ptcl - 20:
        rec_traj[-1] = np.array(rec_traj[-1]).T


rec_surf = []
for temp_idx in mesh.surf_set:
    temp_svec, temp_stheta = mesh.calc_surf_norm(temp_idx, 
                                                  radius=surf_norm_range, 
                                                  imode=surf_norm_mode)
    rec_surf.append([temp_idx, temp_svec])

colMap = copy(cm.Accent)
colMap.set_under(color='white')

def plot_traj(ax, traj):
    # for i in range(10):
    #     ax.plot(traj[num_ptcl - i - 1][0, :], traj[num_ptcl - i - 1][1, :],
    #             marker='o', markersize=0.3, linestyle='-', linewidth=0.1)
    for temp_ptcl in traj:
        ax.plot(temp_ptcl[0, :], temp_ptcl[1, :],
                marker='o', markersize=0.3, linestyle='-', linewidth=0.1)

def plot_surf_norm(ax, posn, svec):
    ax.quiver(posn[0], posn[1],
              svec[0], svec[1], 
              # scale=50, units='xy',
              # headwidth=1, headlength=1, lw=0.01, edgecolors='k',
              width=0.001)


def plot_mesh(mat, ith):
    fig, axes = plt.subplots(1, 2, figsize=(16, 8),
                              constrained_layout=True)
    
    ax = axes[0]
    ax.contourf(mesh.x, mesh.z, mat, cmap=colMap, vmin=0.2, extend='both')
    ax.set_xlim(0.0, mesh.width)
    ax.set_ylim(0.0, mesh.height)
    plot_traj(ax, rec_traj)
    
    ax = axes[1]
    ax.scatter(mesh.x, mesh.z, c=mat, s=1, cmap=colMap, vmin=0.2)
    ax.set_xlim(0.0, mesh.width)
    ax.set_ylim(0.0, mesh.height)
    plot_traj(ax, rec_traj)
    for item in rec_surf:
        temp_idx, temp_svec = item
        temp_posn = np.array([mesh.x[temp_idx], mesh.z[temp_idx]])
        plot_surf_norm(ax, temp_posn, temp_svec)
    
    fig.savefig('mat_%d.png' % ith, dpi=300)

for ith, mat in enumerate(rec_mesh):
    plot_mesh(mat, ith+1)
