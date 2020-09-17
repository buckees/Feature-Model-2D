# -*- coding: utf-8 -*-
"""
Feature Model 2D
Particle Movement
"""

import numpy as np
import math

def find_intersect_node(line, node_surf):
    node_itsct = []
    temp, nptcl = line.shape
    temp, nsurf = node_surf.shape
    for i in range(nptcl):
        a, b, c = line[:, i]
        dist = []
        for j in range(nsurf):
            x0, y0 = node_surf[:, j]
            dist.append(abs(a*x0 + b*y0 + c)/math.sqrt(a**2 + b**2))
        idx = np.argmin(dist)
        node_itsct.append(node_surf[:, idx])
    return np.array(node_itsct).T
            

