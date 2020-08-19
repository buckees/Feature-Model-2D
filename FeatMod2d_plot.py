# -*- coding: utf-8 -*-
"""
Feature Model 2D
Result Plot
"""


def plot_traj(ax, point1, point2):
    ax.plot([point1[0], point2[0]], [point1[1], point2[1]], 'ro-')
