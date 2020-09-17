# -*- coding: utf-8 -*-
"""
Feature Model 2D
Geometry file
"""
import numpy as np
import math

domain_x, domain_z = 100.0, 500.0 # nm
resolution = 1.0 # nm

def create_mat():
    materials = []
    m = [('SiO2', 1), 'rect', (0.0, 0.0, 100.0, 50.0)]
    materials.append(m)
    m = [('Si', 2),   'rect', (0.0, 50.0, 100.0, 300.0)]
    materials.append(m)
    m = [('PR', 3),   'rect', (0.0, 350.0, 30.0, 100.0)]
    materials.append(m)
    m = [('PR', 3),   'rect', (70.0, 350.0, 30.0, 100.0)]
    materials.append(m)
    return materials

materials = create_mat()


def rect_conv(coord):
    ncoord = [math.ceil(item/resolution) for item in coord]
    return [ncoord[0], ncoord[0]+ncoord[2], ncoord[1], ncoord[1]+ncoord[3]]

def find_surface(mesh):
    surf = []
    nx, nz = mesh.shape
    mesh_surf = np.zeros_like(mesh)
    for i in range(nx):
        for j in range(1, nz-1):
            if mesh[i,j]:
                temp = mesh[(i-1) % nx,j]*mesh[(i+1) % nx,j] \
                                *mesh[i,(j-1)% nz]*mesh[i,(j+1) % nz]
                if not temp:
                    mesh_surf[i,j] = 1
                    surf.append((i,j))
            else:
                temp = mesh[(i-1) % nx,j]+mesh[(i+1) % nx,j] \
                                +mesh[i,(j-1)% nz]+mesh[i,(j+1) % nz]
                if temp:
                    mesh_surf[i,j] = 2
    
    for i in [0, nx-1]:
        for j in range(1, nz-1):
            if not mesh[i,j]:
                mesh_surf[i,j] = -1
                surf.append((i,j))
    
    return np.array(surf).T
    