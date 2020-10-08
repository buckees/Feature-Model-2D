"""Feature Model 2D. Mesh file."""

import numpy as np
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

        for material in materials:
            mater = material[0]
            print(mater)
            itype = material[1]
            if itype == 'rect':
                coord = material[2]
                ncoord = rect_conv(coord, self.res_x, self.res_z)
                self.mat[ncoord[2]:ncoord[3],
                         ncoord[0]:ncoord[1]] = int(mater[1])

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
                         * self.mat[(j+1) % self.nz, i]
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
        """Caculate surface normal."""
        temp_mat_surf = self.mat_surf[idx[0]-radius:idx[0]+radius+1,
                                      idx[1]-radius:idx[1]+radius+1]
        temp_x = self.x[idx[0]-radius:idx[0]+radius+1,
                        idx[1]-radius:idx[1]+radius+1]
        temp_z = self.z[idx[0]-radius:idx[0]+radius+1,
                        idx[1]-radius:idx[1]+radius+1]
        print(temp_mat_surf, '\n', temp_x, '\n', temp_z, '\n')

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
        print(min_norm)
        theta = min_norm.x[0] + np.pi/2
        surf_norm = (np.cos(theta), np.sin(theta))
        temp_posn = np.array([self.x[idx], self.z[idx]])
        temp_posn += self.res*surf_norm
        # make sure the surf_norm points out of material
        temp_mat, temp_idx = self.hit_check(temp_posn)
        if temp_mat:
            theta += np.pi
            surf_norm = (np.cos(theta), np.sin(theta))

        return surf_norm


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
    surf_norm = mesh.calc_surf_norm((174, 25))
    print(surf_norm)
