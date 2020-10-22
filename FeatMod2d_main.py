"""Feature Model 2D. Main program."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy

from FeatMod2d_ops import width, height, res_x, res_z, num_ptcl, \
                          threshold, max_rflct
from FeatMod2d_mesh import MESHGRID
# from FeatMod2d_ptcl import PARTICLE
from Species import Arp
from FeatMod2d_rflct import REFLECT

# create mesh
mesh = MESHGRID(width, height, res_x, res_z)
print(mesh)
mesh.mat_input()
mesh.find_surf()
# mesh.plot()

rec_traj = [[] for i in range(num_ptcl)]
rec_surf = []

delta_L = min(res_x, res_z)*0.5
# Arp = PARTICLE('Ar+', 'Ion',  32.0,     1)

Arp_rflct = REFLECT()

for k in range(num_ptcl):
    Arp.dead = 0
    Arp.init_posn(width, height)
    # record initial position
    rec_traj[k].append(Arp.posn.copy())
    Arp.init_uvec('Normal')
#    Arp.init_plot()

    num_rflct = 0
#    while imove_ptcl == 1 and num_rflct < 5:
    for i in range(10000):
        # advance the ptcl by delta_L
        Arp.move_ptcl(delta_L)
        # periodic b.c. at left and right bdry
        Arp.bdry_check(mesh.width, mesh.height, 'periodic')
        # rec_traj ptcl trajectory
        # rec_traj[k].append(Arp.posn.copy())
        if Arp.dead:
            # record ptcl posn when dead
            rec_traj[k].append(Arp.posn.copy())
            break
        hit_mat, hit_idx = mesh.hit_check(Arp.posn)
        if hit_mat:
            # record the hit point
            rec_traj[k].append(Arp.posn.copy())
            # at this position, th ptcl hits a mat
            # decide wehter a reflection or reaction
            mat_name = mesh.mater[hit_mat]
            rnd = np.random.uniform(0.0, 1.0)
            rflct = REFLECT(Arp.name, mat_name, 1.0)
            prob = rflct.calc_prob()
            if rnd < prob:
                if num_rflct > max_rflct:
                    Arp.dead = 1
                    # record ptcl posn when dead
                    rec_traj[k].append(Arp.posn.copy())
                    break
                # call reflection
                Arp_rflct.svec, Arp_rflct.stheta = \
                    mesh.calc_surf_norm(hit_idx, radius=1, imode='Sum Vector')
                rec_surf.append([hit_idx, Arp_rflct.svec.copy()])
                # use only specular reflection
                Arp.uvec = Arp_rflct.spec_rflct(Arp.uvec)
                # Arp.uvec = Arp_rflct.rflct(Arp.uvec)
                # update ptcl position as the hit cell center
                # Arp.posn = np.array([mesh.x[hit_idx], mesh.z[hit_idx]])
                num_rflct += 1
            else:
                # now ireact = 1
                mesh.update_mat(hit_idx, threshold)
                mesh.find_float_cell()
                mesh.find_surf()
                Arp.dead = 1
        # check if the ptcl is dead
        if Arp.dead:
            break

    rec_traj[k] = np.array(rec_traj[k]).T

for temp_idx in mesh.surf:
    temp_svec, temp_stheta = mesh.calc_surf_norm(temp_idx)
    rec_surf.append([temp_idx, temp_svec])

colMap = copy.copy(cm.Accent)
colMap.set_under(color='white')

def plot_traj(ax, traj):
    for i in range(10):
        ax.plot(traj[num_ptcl - i - 1][0, :], traj[num_ptcl - i - 1][1, :],
                marker='o', markersize=0.3, linestyle='-', linewidth=0.1)

def plot_surf_norm(ax, posn, svec):
    ax.quiver(posn[0], posn[1],
              svec[0], svec[1])


fig, axes = plt.subplots(1, 2, figsize=(4, 8),
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
for item in rec_surf:
    temp_idx, temp_svec = item
    temp_posn = np.array([mesh.x[temp_idx], mesh.z[temp_idx]])
    plot_surf_norm(ax, temp_posn, temp_svec)

plt.show()
fig.savefig('etching_demo.png', dpi=600)
