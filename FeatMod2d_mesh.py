"""Feature Model 2D. Mesh file."""

import numpy as np
from math import cos, sin
import copy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import minimize, least_squares


class MESHGRID(object):
    """Mesh object."""

    def __init__(self, width, height, res_x, res_z):
        self.width = width  # domain width in x direction
        self.height = height  # domain height in z direction
        self.res_x = res_x  # resolution in x direction
        self.res_z = res_z  # resolution in z direction
        # array of resolution
        self.res = np.array([self.res_x, self.res_z])
#        <--------------------- width --------------------->
#        0.0 ---- 1.0 ---- 2.0 ---- ... ---- 99.0 ---- 100.0
#        --- cell --- cell --- cell ... cell ---- cell -----
#        -- center - center - center . center -- center ----
#        --- 0.5 ---- 1.5 ---- 2.5  ... 98.5 ---- 99.5 -----
#        --(0, nj)--(1, nj)--(2, nj)...(98, nj)--(99, nj)---
        self.nx = int(np.ceil(width/res_x))  # num of cells in x
        self.nz = int(np.ceil(height/res_z))  # num of cels iin z
        # init x and z coordinates
        tempx = np.linspace(0.0, self.width, self.nx)
        tempz = np.linspace(0.0, self.height, self.nz)
        self.x, self.z = np.meshgrid(tempx, tempz)
        # mesh materials is assigned to self.mat matrix
        # note that the shape of self.mat is (nz, nx)
        self.mat = np.zeros_like(self.x).astype(int)
        # construct a mat-like matrix for surface
        # mat_surf = 1 if a surface node; 0 if not.
        self.mat_surf = np.zeros_like(self.mat).astype(int)
        self.surf = np.array([])  # surface node
        self.mater = []  # materials name <--> materails No.

    def __str__(self):
        """Print out mesh information."""
        return """
               This mesh with domain of %.1f nm (width) x %.1f nm (height)
               and resolution in (x, z) = (%.2f nm, %.2f nm)
               and number of cells in (x, z) = (%d, %d)
               """ \
                % (self.width, self.height,
                   self.res_x, self.res_z,
                   self.nx, self.nz)

    def mat_input(self):
        """Assign input materials to the mesh."""
        self.mater = ['Vac', 'SiO2', 'Si', 'PR']
        materials = []
        m = [('SiO2', 1), 'rect', (0.0, 0.0, 100.0, 50.0)]
        materials.append(m)
        m = [('Si', 2),   'rect', (0.0, 50.0, 100.0, 300.0)]
        materials.append(m)
        m = [('PR', 3),   'rect', (0.0, 350.0, 30.0, 100.0)]
        materials.append(m)
        m = [('PR', 3),   'rect', (70.0, 350.0, 30.0, 100.0)]
        materials.append(m)
        # m = [('Plasma', 0),   'circ', (50.0, 350.0, 30.0)]
        # materials.append(m)

        for material in materials:
            mater = material[0]
            print(mater)
            itype = material[1]
            if itype == 'rect':
                coord = material[2]
                ncoord = rect_conv(coord, self.res_x, self.res_z)
                self.mat[ncoord[2]:ncoord[3],
                         ncoord[0]:ncoord[1]] = int(mater[1])
            if itype == 'circ':
                circ_x, circ_z, circ_r = material[2]
                for j in range(1, self.nz-1):
                    for i in range(self.nx):
                        tempx = self.x[j, i]
                        tempz = self.z[j, i]
                        tempd = (tempx - circ_x)**2 + (tempz - circ_z)**2
                        if tempd < circ_r**2:
                            self.mat[j, i] = int(mater[1])

    def find_surf(self):
        """Search for the surface nodes."""
        surf = []
        # search surface within materials
        for j in range(1, self.nz-1):
            for i in range(self.nx):
                # if mat[i,j] is not 0
                if self.mat[j, i]:
                    # temp = multiply all neighbours
                    # temp = 0 means one of neighbours = 0
                    temp = self.mat[j, (i-1) % self.nx] \
                         * self.mat[j, (i+1) % self.nx] \
                         * self.mat[(j-1) % self.nz, i] \
                         * self.mat[(j+1) % self.nz, i] \
                         * self.mat[(j-1) % self.nz, (i-1) % self.nx] \
                         * self.mat[(j-1) % self.nz, (i+1) % self.nx] \
                         * self.mat[(j+1) % self.nz, (i-1) % self.nx] \
                         * self.mat[(j+1) % self.nz, (i+1) % self.nx]
                    if not temp:
                        surf.append((j, i))
                        self.mat_surf[j, i] = 1

        # search for left and right boundaries
#        for j in range(1, self.nz-1):
#            for i in [0, self.nx-1]:
#                if not self.mat[j,i]:
#                    surf.append((j,i))

        self.surf = np.array(surf).T

    def plot(self):
        """Plot mesh and surface."""
        colMap = copy.copy(cm.get_cmap("Accent"))
        colMap.set_under(color='white')

        fig, axes = plt.subplots(1, 2, figsize=(4, 8),
                                 constrained_layout=True)
        ax = axes[0]
        ax.contourf(self.x, self.z, self.mat,
                    cmap=colMap, vmin=0.2, extend='both')
        ax = axes[1]
        ax.contourf(self.x, self.z, self.mat_surf,
                    cmap=colMap, vmin=0.2, extend='both')
        # axes[1].plot(self.surf[1, :]*self.res_x,
        #              self.surf[0, :]*self.res_z,
        #              'o', )
        # axes[1].set_xlim(0, self.width)
        # axes[1].set_ylim(0, self.height)
        plt.show(fig)

    def hit_check(self, posn):
        """
        Check wether a particle hits a material.

        All posn needs to be rounded to 0.5*nx (cell center).
        posn needs to be shifted by half res_x.
        """
        idx = np.rint((posn - self.res*0.5) / self.res).astype(int)
        # reverse idx in order to accomplish self.mat order
        # print('hit = ', posn, idx)
        idx = np.flipud(idx)
        # convert idx to index format
        idx = tuple(idx)
        return self.mat[idx], idx

    def update_mat(self, idx, threshold):
        """
        Update materials due to etching.

        threshold: etching probability
        mat == 2 equivalent to mat == 'Si'
        """
        if self.mat[idx] == 2:
            rnd = np.random.uniform(0.0, 1.0)
            if rnd < threshold:
                self.mat[idx] = 0

    def calc_surf_norm(self, idx, radius=2):
        """
        Caculate surface normal.

        Inputs: index where particle hit
        Calc searching sub-domain in which the surface fits
        Calc cost function for surface fitting
        Calc the minimun cost function
        Calc the surface direction
        Calc the surface normal direction, rotate 90 or 270 degrees
        Calc surface normal vector
        Output: surface normal vector and vector angle
        """
        # Calc the sub-domain boundary
        bottom = idx[0]-radius
        top = idx[0]+radius+1
        left = idx[1]-radius
        right = idx[1]+radius+1
        # cut sub-domain if it is out of main-domain
        if bottom < 0:
            bottom = 0
        if top > self.nz-1:
            top = self.nz-1
        if left < 0:
            left = 0
        if right > self.nx-1:
            right = self.nx-1
        # Construct the sub-domain
        temp_mat_surf = self.mat_surf[bottom:top, left:right]
        temp_x = self.x[bottom:top, left:right]
        temp_z = self.z[bottom:top, left:right]

        def cost_func_surf_norm(theta):
            """Construct the cost func for surface fitting."""
            A, B = -np.sin(theta), np.cos(theta)
            C = A*self.x[idx] + B*self.z[idx]
            Q = A*temp_x + B*temp_z - C
            Q = np.multiply(Q, temp_mat_surf)
            Q2 = np.power(Q, 2)
            Qsum = Q2.sum()
            return Qsum

        min_norm = minimize(cost_func_surf_norm, np.pi/4)
        # obtain the surface normal
        theta = min_norm.x[0] + np.pi/2
        surf_norm = np.array([cos(theta), sin(theta)])
        temp_posn = np.array([self.x[idx], self.z[idx]])
        temp_posn += np.sqrt(2)/2*radius*self.res*surf_norm
        # Check boundaries
        temp_posn[0] = np.clip(temp_posn[0], 0.0, self.width-self.res_x*1e-3)
        temp_posn[1] = np.clip(temp_posn[1], 0.0, self.height-self.res_z*1e-3)
        # make sure the surf_norm points out of material
        temp_mat, temp_idx = self.hit_check(temp_posn)
        if temp_mat:
            theta += np.pi
            surf_norm = np.array([cos(theta), sin(theta)])

        return surf_norm, theta

    def find_float_cell(self):
        """
        Search for the cells which are not connected.
        
        Scan each cell, check its 4 neighours, top, bottom, left and right.
        If they are all empty, the cell is identified as a floating cell, 
        which will be dropped to bottom.
        """
        for j in range(1, self.nz-1):
            for i in range(self.nx):
                # if mat[i,j] is not 0
                if self.mat[j, i]:
                    # check its 4 neighbours
                    bottom = self.mat[j-1, i]
                    top = self.mat[j+1, i]
                    if i == 0:
                        left = self.mat[j, self.nx-1]
                    else:
                        left = self.mat[j, i-1]
                    if i == self.nx-1:
                        right = self.mat[j, 0]
                    else:
                        right = self.mat[j, i+1]
                    if (bottom + top + left + right) == 0:
                        self.drop_cell((j, i))

    def drop_cell(self, idx):
        """
        Drop the cells at idx.
        
        Drop the cell by 1 cell down until it hit bottom materials.
        """
        idx_j, idx_i = idx
        # remove the cell at idx
        temp_mat = self.mat[idx]
        self.mat[idx] = 0
        bottom = self.mat[idx_j-1, idx_i]
        while bottom == 0:
            idx_j -= 1
            bottom = self.mat[idx_j-1, idx_i]
        self.mat[idx_j, idx_i] = temp_mat

def rect_conv(coord, res_x, res_z):
    """Convert rectangular coordinates to index."""
    ncoord = [np.ceil(coord[0]/res_x), np.ceil(coord[1]/res_z),
              np.ceil(coord[2]/res_x), np.ceil(coord[3]/res_z)]
    ncoord = [int(temp) for temp in ncoord]
    return [ncoord[0], ncoord[0]+ncoord[2], ncoord[1], ncoord[1]+ncoord[3]]


if __name__ == '__main__':
    from FeatMod2d_ops import width, height, res_x, res_z
    mesh = MESHGRID(width, height, res_x, res_z)
    mesh.mat_input()
    mesh.find_surf()
    mesh.plot()
    temp_idx = tuple(mesh.surf[:, 1])
    print('idx=', temp_idx)
    surf_norm, theta = mesh.calc_surf_norm(temp_idx)
    print(surf_norm)
    temp_posn = temp_idx*mesh.res

    colMap = copy.copy(cm.get_cmap("Accent"))
    colMap.set_under(color='white')

    fig, axes = plt.subplots(1, 2, figsize=(4, 8),
                             constrained_layout=True)
    ax = axes[0]
    ax.contourf(mesh.x, mesh.z, mesh.mat,
                cmap=colMap, vmin=0.2, extend='both')
    ax.quiver(temp_posn[1], temp_posn[0], surf_norm[0], surf_norm[1],
              scale=5)
    ax = axes[1]
    ax.contourf(mesh.x, mesh.z, mesh.mat_surf,
                cmap=colMap, vmin=0.2, extend='both')
    plt.show(fig)
