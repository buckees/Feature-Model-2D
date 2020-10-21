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

delta_L = min(res_x, res_z)
# Arp = PARTICLE('Ar+', 'Ion',  32.0,     1)

Arp_rflct = REFLECT()

for k in range(num_ptcl):
    Arp.dead = 0
    Arp.init_posn(width, height)
    Arp.init_uvec('Normal')
#    Arp.init_plot()
    rec_traj[k].append(Arp.posn.copy())

    num_rflct = 0
#    while imove_ptcl == 1 and num_rflct < 5:
    for i in range(10000):
        # advance the ptcl by delta_L
        Arp.move_ptcl(delta_L)
        # periodic b.c. at left and right bdry
        Arp.bdry_check(mesh.width, mesh.height, 'periodic')
        # rec_traj ptcl trajectory
        rec_traj[k].append(Arp.posn.copy())
        if Arp.dead:
            break
        hit_mat, hit_idx = mesh.hit_check(Arp.posn)
        if hit_mat:
            # at this position, th ptcl hits a mat
            # decide wehter a reflection or reaction
            rnd = np.random.uniform(0.0, 1.0)
            mat_name = mesh.mater[hit_mat]
            rflct = REFLECT(Arp.name, mat_name, 1.0)
            prob = rflct.calc_prob()
            if rnd < prob:
                if num_rflct > max_rflct:
                    Arp.dead = 1
                    break
                # call reflection
                Arp_rflct.svec, Arp_rflct.stheta = mesh.calc_surf_norm(hit_idx)
                rec_surf.append([hit_idx, Arp_rflct.svec.copy()])
                Arp.uvec = Arp_rflct.rflct(Arp.uvec)
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


colMap = copy.copy(cm.Accent)
colMap.set_under(color='white')

def plot_traj(ax, traj):
    for i in range(10):
        ax.plot(traj[num_ptcl - i - 1][0, :], traj[num_ptcl - i - 1][1, :],
                marker='o', markersize=1, linestyle='-')

def plot_surf_norm(ax, posn, svec):
    ax.quiver(posn[1], posn[0],
              svec[0], svec[1])

def idx2posn(idx):
    """
    Convert index to position.

    idx: tuple (j, i), where j is vertical index and i is horizontal index.
    """
    idx = np.array(idx)
    posn = np.multiply(idx, mesh.res) + mesh.res*0.5
    return posn

print(rec_surf)

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
    temp_posn = idx2posn(temp_idx)
    plot_surf_norm(ax, temp_posn, temp_svec)

plt.show()
fig.savefig('etching_demo.png', dpi=600)
