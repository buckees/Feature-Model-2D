# -*- coding: utf-8 -*-
"""
Feature Model 2D
Result Plot
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from FeatMod2d_geo import domain_x, domain_z
from FeatMod2d_ops import num_ptcl

colMap = cm.Accent
colMap.set_under(color='white')

def plot_mesh(mesh, surf):            
    fig, axes = plt.subplots(1,2, figsize=(4,8),
                               constrained_layout=True)
    axes[0].contourf(mesh.T, cmap = colMap, vmin = 0.2, extend='both')
    axes[1].plot(surf[0, :], surf[1, :], 'o', )
    axes[1].set_xlim(0, domain_x)
    axes[1].set_ylim(0, domain_z)
    plt.show(fig)

def plot_init():
    fig, axes = plt.subplots(1,2, figsize=(8,3),
                           constrained_layout=True)
    axes[0].plot(ux, uz, 'o') 
    axes[1].quiver(0.0, 0.0, ux, uz)
    plt.show(fig2)

def plot_itsct(mesh, posn_init, node_traj, node_itsct):
    fig, axes = plt.subplots(1,2, figsize=(4,8),
                               constrained_layout=True)
    axes[0].contourf(mesh.T, cmap = colMap, vmin = 0.2, extend='both')
    axes[0].plot(node_itsct[0, :], node_itsct[1, :], 'ro')    
    for i in range(num_ptcl):
        point1 = posn_init[:, i]
        point2 = node_traj[:, i]
        plot_traj(axes[0], point1, point2, 'green')
    axes[0].set_xlim(0, domain_x)
    axes[0].set_ylim(0, domain_z)
    axes[1].contourf(mesh.T, cmap = colMap, vmin = 0.2, extend='both')
    for i in range(num_ptcl):
        point1 = posn_init[:, i]
        point2 = node_itsct[:, i]
        plot_traj(axes[1], point1, point2, 'red')
    axes[1].set_xlim(0, domain_x)
    axes[1].set_ylim(0, domain_z)
    plt.show(fig)

def plot_traj(ax, point1, point2, color):
    ax.plot([point1[0], point2[0]], [point1[1], point2[1]], 
            color=color, marker='o', linestyle='-')
