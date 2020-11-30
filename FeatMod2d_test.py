"""Feature Model 2D. Main program."""

import os
import glob
for i in glob.glob("*.png"):
    os.remove(i)

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from copy import copy, deepcopy
import numpy as np
from random import choices


from FeatMod2d_ops import (width, height, res_x, res_z, num_ptcl, ibc, 
                          threshold, max_rflct, idstrb, step_fac, max_step,
                          num_plot, surf_norm_range, surf_norm_mode)
from FeatMod2d_mesh import MESHGRID
from FeatMod2d_ptcl import PARTICLE
from FeatMod2d_readin import (sp_run_list, sp_name, sp_weight, mat_name)
from FeatMod2d_rct import React
from FeatMod2d_mat import Si2d, Si2d_trench_v02
from FeatMod2d_rflct import REFLECT

# create mesh
mesh = MESHGRID(width, height, res_x, res_z)
print(mesh)
mesh.add_mat(Si2d_trench_v02)

# Frame the chem data
react = React(sp_name=sp_name, mat_name=mat_name)
react.readin('Chem.csv')
print(react.df)

delta_L = min(res_x, res_z)*step_fac

# Initialize the PARTICLE() object
# ptcl = PARTICLE()
ptcl_rflct = REFLECT()

# record for diagnostics
rec_traj = []

for k in range(num_ptcl):
    # randomly choose species from the sp_run_list according to its flux weight
    temp_sp, = choices(sp_name, weights=sp_weight, k=1)
    temp_sp_info = sp_run_list[temp_sp]
    ptcl = PARTICLE(**temp_sp_info)
    if (k + 1) % int(num_ptcl/num_plot) == 0:
        print('%d particles are launched!' % (k+1))
        mesh.plot(dpi=300, fname='nptcl=%d.png' % (k+1))
        mesh.plot_surf(surf_norm_range=surf_norm_range, 
                       surf_norm_mode=surf_norm_mode, 
                       dpi=300, fname='surf_nptcl=%d.png' % (k+1))
        # mesh.find_float_cell(idiag=1)

    ptcl.dead = 0
    ptcl.init_posn(width, height)
    if ptcl.ptype == 'Ion':
        idstrb=['Uniform2D', -5.0, 5.0]
    elif ptcl.ptype == 'Neut':
        idstrb=['Uniform2D', -45.0, 45.0]
    ptcl.init_uvec(idstrb)
    num_rflct = 0
    
    # record the particle trajectory
    rec_traj.append([])
    rec_traj[-1].append(ptcl.posn.copy())

    for i in range(max_step):
        # advance the ptcl by delta_L
        ptcl.move_ptcl(delta_L)
        # periodic b.c. at left and right bdry
        ptcl.bdry_check(mesh.width, mesh.height, ibc)
        # check if the ptcl is dead
        if ptcl.dead:
            rec_traj[-1].append(ptcl.posn.copy())
            break
        hit_mat, hit_idx = mesh.hit_check(ptcl.posn)
        if hit_mat:
            rec_traj[-1].append(ptcl.posn.copy())
            ptcl_rflct.svec, ptcl_rflct.stheta = \
                mesh.calc_surf_norm(hit_idx, radius=surf_norm_range, 
                                    imode=surf_norm_mode)
            
            # determine react or rflct
            hit_sp_name = ptcl.name
            hit_mat_name = mesh.mater[hit_mat]
            hit_df = react.df[
                    (react.df['sp'] == hit_sp_name) 
                    & (react.df['mat'] == hit_mat_name)
                    ]
            
            row, col = hit_df.shape
            
            if row > 1:
                # now more than rflct can occur
                hit_df_idx = list(hit_df.index)
                hit_df_prob = list(hit_df['prob'].tolist())    
                chosen_idx, = choices(hit_df_idx, weights=hit_df_prob, k=1)
                chosen_df = hit_df.loc[chosen_idx]
                if chosen_df['type'] == 'rflct':
                    break
                elif chosen_df['type'] == 'etch':
                    mesh.mat[hit_idx] = 0
                    mesh.update_surf(hit_idx)
                    ptcl.dead = 1
                elif chosen_df['type'] == 'chem':
                    mesh.mat[hit_idx] = 4
                    ptcl.dead = 1
            else:
                if num_rflct > max_rflct:
                    ptcl.dead = 1
                    rec_traj[-1].append(ptcl.posn.copy())
                    break
                # reflect
                ptcl.uvec = ptcl_rflct.rflct(ptcl.uvec)
                # move the ptcl by 10 steps until it gets out of the mat
                for ii in range(10):
                    ptcl.move_ptcl(delta_L)
                    ptcl.bdry_check(mesh.width, mesh.height, ibc)
                    hit_mat, hit_idx = mesh.hit_check(ptcl.posn)
                    if not hit_mat:
                        break
                if hit_mat:
                    ptcl.dead = 1
                    rec_traj[-1].append(ptcl.posn.copy())
                    break
                num_rflct += 1
                
        # check if the ptcl is dead
        if ptcl.dead:
            rec_traj[-1].append(ptcl.posn.copy())
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

np.save('mesh_mat', mesh.mat)
np.save('mesh_surf', mesh.surf)
np.save('mesh_x', mesh.x)
np.save('mesh_z', mesh.z)

np.savetxt('mesh_mat.csv', mesh.mat, delimiter=',')
np.savetxt('mesh_surf.csv', mesh.surf, delimiter=',')
np.savetxt('mesh_x.csv', mesh.x, delimiter=',')
np.savetxt('mesh_z.csv', mesh.z, delimiter=',')
