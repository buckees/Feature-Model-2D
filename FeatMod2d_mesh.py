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
        # init x, z coordinates, which represent the cell center 
        tempx = np.linspace(0.0 + 0.5*self.res_x, 
                            self.width - 0.5*self.res_x, self.nx)
        tempz = np.linspace(0.0 + 0.5*self.res_z, 
                            self.height - 0.5*self.res_z, self.nz)
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
        # m = [('Vac', 0),   'circ', (50.0, 350.0, 30.0)]
        # materials.append(m)

        for material in materials:
            mater = material[0]
            # print(mater)
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
        self.surf = []
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
                        self.surf.append((j, i))
                        self.mat_surf[j, i] = 1

        # search for left and right boundaries
#        for j in range(1, self.nz-1):
#            for i in [0, self.nx-1]:
#                if not self.mat[j,i]:
#                    surf.append((j,i))

        # self.surf = np.array(surf).T

    def find_surf_vac(self):
        """Search for the surface nodes neighbors in vac."""
        self.surf_vac = []
        # search surface within materials
        for j in range(1, self.nz-1):
            for i in range(self.nx):
                # if mat[i,j] is 0
                if not self.mat[j, i]:
                    # temp = sum all neighbours
                    # temp != 0 means one of neighbours is surf node
                    temp = self.mat[j, (i-1) % self.nx]
                    temp += self.mat[j, (i+1) % self.nx]
                    temp += self.mat[(j-1) % self.nz, i]
                    temp += self.mat[(j+1) % self.nz, i]
                    temp += self.mat[(j-1) % self.nz, (i-1) % self.nx]
                    temp += self.mat[(j-1) % self.nz, (i+1) % self.nx]
                    temp += self.mat[(j+1) % self.nz, (i-1) % self.nx]
                    temp += self.mat[(j+1) % self.nz, (i+1) % self.nx]
                    if temp:
                        self.surf_vac.append((j, i))
                        self.mat_surf[j, i] = -1


    def plot(self):
        """Plot mesh and surface."""
        colMap = copy.copy(cm.get_cmap("Accent"))
        colMap.set_under(color='white')

        fig, axes = plt.subplots(1, 2, figsize=(4, 8), dpi=600,
                                 constrained_layout=True)
        ax = axes[0]
        ax.scatter(self.x, self.z, c=self.mat, s=1, cmap=colMap, vmin=0.2)
        ax = axes[1]
        ax.scatter(self.x, self.z, c=self.mat_surf, s=1)
        
    def hit_check(self, posn):
        """
        Check wether a particle hits a material.

        All posn needs to be rounded to 0.5*nx (cell center).
        posn needs to be shifted by half res_x.
        """
        idx = np.rint((posn - self.res*0.5) / self.res).astype(int)
        # reverse idx in order to accomplish self.mat order
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

    def calc_surf_norm(self, idx, radius=1, imode="Fit Plane", bc='Periodic'):
        """
        Caculate surface normal.

        idx: index where particle hit
        radius: the sub-domain where the surf norm is calc
        imode = 1: fitting plane to the surf sites within the sub-domain
                2: sum the vector of hit site to vac sites 
        
        Search for the sub-domain in which the surface fits
        imode = 1
        Calc cost function for surface fitting
        Calc the minimun cost function
        Calc the surface direction
        Calc the surface normal direction, rotate 90 or 270 degrees
        Calc surface normal vector
        
        imode = 2
        Calc the vec from hit surf node to surf_vac nodes
        Sum all the vec (weight can be specified)
        Normal the vec
        
        Output: surface normal vector and vector angle
        """
        # Check input
        if imode in ['Fit Plane', 'Sum Vector']:
            pass
        else:
            return print('Error')
        # Create the sub-domain boundary
        bottom = idx[0]-radius
        top = idx[0]+radius+1
        left = idx[1]-radius
        right = idx[1]+radius+1
        # cut sub-domain if it is out of main-domain vertically
        if bottom < 0:
            bottom = 0
        if top > self.nz-1:
            top = self.nz-1
        # if left < 0:
        #     left = 0
        # if right > self.nx-1:
        #     right = self.nx-1
        # Construct the sub-domain, periodic b.c.
        if left < 0:
            sub_mat_surf = np.concatenate((self.mat_surf[bottom:top, left:], 
                           self.mat_surf[bottom:top, 0:right]), axis=1)
            sub_x = np.concatenate((self.x[bottom:top, left:] - self.width, 
                                    self.x[bottom:top, 0:right]), 
                                   axis=1)
            sub_z = np.concatenate((self.z[bottom:top, left:], 
                                    self.z[bottom:top, 0:right]), axis=1)
        elif right > self.nx-1:
            right = right % self.nx
            sub_mat_surf = np.concatenate((self.mat_surf[bottom:top, left:], 
                           self.mat_surf[bottom:top, 0:right]), axis=1)
            sub_x = np.concatenate((self.x[bottom:top, left:], 
                                    self.x[bottom:top, 0:right] + self.width), 
                                   axis=1)
            sub_z = np.concatenate((self.z[bottom:top, left:], 
                                    self.z[bottom:top, 0:right]), axis=1)
        else:
            sub_mat_surf = self.mat_surf[bottom:top, left:right]
            sub_x = self.x[bottom:top, left:right]
            sub_z = self.z[bottom:top, left:right]
        # print(sub_mat_surf)

        if imode == "Fit Plane":
            # sub_mat_surf consists of 1(surf) and -1(surf_vac)
            # surf_vac is not used when imode == 1, zero out -1
            temp_sub_mat_surf = np.where(sub_mat_surf == -1, 0, sub_mat_surf)
            def cost_func_surf_norm(theta):
                """Construct the cost func for surface fitting."""
                A, B = -np.sin(theta), np.cos(theta)
                C = A*self.x[idx] + B*self.z[idx]
                Q = A*sub_x + B*sub_z - C
                Q = np.multiply(Q, temp_sub_mat_surf)
                Qsum = -abs(Q.sum())
                Q2 = np.power(Q, 2)
                Qsum += Q2.sum()
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
        
        elif imode == "Sum Vector":
            # sub_mat_surf consists of 1(surf) and -1(surf_vac)
            # surf_vac is not used when imode == 2, zero out 1
            temp_sub_mat_surf = np.where(sub_mat_surf == 1, 0, sub_mat_surf)
            # print(temp_sub_mat_surf)
            temp_vecx = np.multiply(self.x[idx] - sub_x, temp_sub_mat_surf)
            temp_vecz = np.multiply(self.z[idx] - sub_z, temp_sub_mat_surf)
            # Calc the vector norm**2
            temp_vec_norm = np.power(temp_vecx, 2) + np.power(temp_vecz, 2)
            # The far from the idx, the smaller the weight. weight = 1/r
            temp_vecx = np.divide(temp_vecx, temp_vec_norm, 
                                  out=np.zeros_like(temp_vec_norm), 
                                  where=temp_vec_norm!=0)
            temp_vecz = np.divide(temp_vecz, temp_vec_norm,
                                  out=np.zeros_like(temp_vec_norm), 
                                  where=temp_vec_norm!=0)
            # print(temp_vecz)
            temp_vecx = temp_vecx.sum()
            temp_vecz = temp_vecz.sum()
            surf_norm = np.array([temp_vecx, temp_vecz])
            temp_norm = np.linalg.norm(surf_norm)
            if temp_norm:
                surf_norm = surf_norm/np.linalg.norm(surf_norm)
                theta = np.arccos(surf_norm)
            else:
                theta = np.random.uniform(-np.pi, np.pi)
                surf_norm = np.array([sin(theta), -cos(theta)])
        
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
    mesh.find_surf_vac()
    mesh.plot()
    
    colMap = copy.copy(cm.get_cmap("Accent"))
    colMap.set_under(color='white')

    rec_surf = []
    for temp_idx in mesh.surf:
        # temp_idx = (224, 13)
        temp_svec, temp_stheta = mesh.calc_surf_norm(temp_idx, 
                                                     radius=3,
                                                     imode='Fit Plane')
        rec_surf.append([temp_idx, temp_svec])

    def plot_surf_norm(ax, posn, svec):
        ax.quiver(posn[0], posn[1],
                  svec[0], svec[1])


    fig, axes = plt.subplots(1, 2, figsize=(4, 8),
                             constrained_layout=True)
    
    ax = axes[0]
    ax.contourf(mesh.x, mesh.z, mesh.mat, cmap=colMap, vmin=0.2, extend='both')
    ax.set_xlim(0.0, mesh.width)
    ax.set_ylim(0.0, mesh.height)
    
    ax = axes[1]
    ax.scatter(mesh.x, mesh.z, c=mesh.mat, s=1, cmap=colMap, vmin=0.2)
    ax.set_xlim(0.0, mesh.width)
    ax.set_ylim(0.0, mesh.height)

    for item in rec_surf:
        temp_idx, temp_svec = item
        temp_posn = np.array([mesh.x[temp_idx], mesh.z[temp_idx]])
        plot_surf_norm(ax, temp_posn, temp_svec)
    
    plt.show()
    fig.savefig('init.png', dpi=600)
