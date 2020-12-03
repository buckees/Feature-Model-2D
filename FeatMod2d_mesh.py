"""Feature Model 2D. Mesh file."""

import numpy as np
from math import cos, sin
import copy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
colMap = copy.copy(cm.get_cmap("Accent"))
colMap.set_under(color='white')

from scipy.optimize import minimize, least_squares

from FeatMod2d_misc import isInsideTrgl


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
        # surf = 1 if a surface node; 0 if not.
        self.surf = np.zeros_like(self.mat).astype(int)
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

    def add_mat(self, materials):
        """Assign input materials to the mesh."""
        self.mater = ['Vac_', 'PR_', 'Si_', 'SiO2_', 'SiCl_']

        for material in materials:
            mater = material[0]
            itype = material[1]
            
            # check itype error
            if itype in ['rect', 'circ', 'trgl']:
                pass
            else:
                return print('itype = %s is not found' % itype)
            
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
                
            if itype == 'trgl':
                p1_x, p1_y, p2_x, p2_y, p3_x, p3_y = material[2]
                for j in range(1, self.nz-1):
                    for i in range(self.nx):
                        tempx = self.x[j, i]
                        tempz = self.z[j, i]
                        if isInsideTrgl(p1_x, p1_y, p2_x, p2_y, p3_x, p3_y,
                                        tempx, tempz):
                            self.mat[j, i] = int(mater[1])
        
        # find the surf nodes and surf_vac nodes
        self._find_surf()
        self._find_surf_set()
    
    def restart_mat(self, fname='restart_mat.npy'):
        """Read in mat from restart file."""
        self.mat = np.load(fname)
        # find the surf nodes and surf_vac nodes
        self._find_surf()
        self._find_surf_set()

    def _check_surf(self, idx):
        j, i = idx
        # if self.mat[idx] > 0, it could be a surf node
        if self.mat[idx]:
            # find surf nodes
            # tempa = multiply all neighbours
            # tempa = 0 means one of neighbours = 0
            tempa = self.mat[j, (i-1) % self.nx]
            tempa *= self.mat[j, (i+1) % self.nx]
            tempa *= self.mat[(j-1) % self.nz, i]
            tempa *= self.mat[(j+1) % self.nz, i]
            tempa *= self.mat[(j-1) % self.nz, (i-1) % self.nx]
            tempa *= self.mat[(j-1) % self.nz, (i+1) % self.nx]
            tempa *= self.mat[(j+1) % self.nz, (i-1) % self.nx]
            tempa *= self.mat[(j+1) % self.nz, (i+1) % self.nx]
            if not tempa:
                self.surf[idx] = 1
        # if self.mat[idx] = 0, it could be a surf_vac node
        else:    
            # find surf nodes in vac
            # tempb = sum all neighbours
            # tempb != 0 means one of neighbours is surf node
            tempb = self.mat[j, (i-1) % self.nx]
            tempb += self.mat[j, (i+1) % self.nx]
            tempb += self.mat[(j-1) % self.nz, i]
            tempb += self.mat[(j+1) % self.nz, i]
            tempb += self.mat[(j-1) % self.nz, (i-1) % self.nx]
            tempb += self.mat[(j-1) % self.nz, (i+1) % self.nx]
            tempb += self.mat[(j+1) % self.nz, (i-1) % self.nx]
            tempb += self.mat[(j+1) % self.nz, (i+1) % self.nx]
            if tempb:
                self.surf[idx] = -1

    def _find_surf(self):
        """Search for the surface nodes."""
        self.surf = np.zeros_like(self.mat).astype(int)
        # search surface within materials
        for j in range(1, self.nz-1):
            for i in range(self.nx):
                self._check_surf((j, i))

    def _find_surf_set(self):
        # construct the surf set
        self.surf_set = set()
        # find the beginning node, search starting from the top left corner
        for j in reversed(range(self.nz)):
            if self.mat[(j, 0)]:
                self.surf_set.add((j, 0))
                _idx_curr = (j, 0)
                break

        # using recursive method for depth search globally
        self._find_next_node(_idx_curr, 
                             lb=0, rb=self.nx-1, tb=self.nz-1, bb=0)
        
    def _find_next_node(self, idx_curr, 
                        lb=0, rb=100, tb=100, bb=0):
        """
        Depth search for surf nodes. 
        
        Find next node, search directions in sequence: left, up, right, down.
        idx_curr: a.u., (j, i) tuple, starting node for depth search.
        lb, rb, tb, bb: a.u., int, left, right, top and bottom bndy.
        """
        _j, _i = idx_curr
        _surf_left, _surf_up, _surf_right, _surf_down = 0, 0, 0, 0
        _idx_left, _idx_up, _idx_right, _idx_down = \
                        (_j, _i-1), (_j+1, _i), (_j, _i+1), (_j-1, _i)
        # search in left
        if (_i-1) >= lb:
            if not (_idx_left in self.surf_set):
                _surf_left = self.surf[_idx_left]
                if _surf_left == 1:
                    self.surf_set.add(_idx_left)
                    self._find_next_node(_idx_left,
                                         lb=lb, rb=rb, tb=tb, bb=bb)
        # search in right    
        if (_i+1) <= rb:
            if not (_idx_right in self.surf_set):
                _surf_right = self.surf[_idx_right]
                if _surf_right == 1:
                    self.surf_set.add(_idx_right)
                    self._find_next_node(_idx_right,
                                         lb=lb, rb=rb, tb=tb, bb=bb)
        # search in up
        if (_j+1) <= tb:
            if not (_idx_up in self.surf_set):
                _surf_up = self.surf[_idx_up]
                if _surf_up == 1:
                    self.surf_set.add(_idx_up)
                    self._find_next_node(_idx_up,
                                         lb=lb, rb=rb, tb=tb, bb=bb)
        # search in down
        if (_j-1) >= bb:
            if not (_idx_down in self.surf_set):
                _surf_down = self.surf[_idx_down]
                if _surf_down == 1:
                    self.surf_set.add(_idx_down)
                    self._find_next_node(_idx_down,
                                         lb=lb, rb=rb, tb=tb, bb=bb)

    
    def update_surf(self, idx, radius=2):
        """
        Search for the surface nodes.
        
        Using find_surf() results in high computational cost.
        Therefore, only the neighbors of the changed node is re-searched.
        """
        _j, _i = idx
        # update the surf nodes within a box of 2*radius+1
        for j in range(_j-radius, _j+radius+1):
            for i in range(_i-radius, _i+radius+1):
                self.surf[j, i] = 0
                self._check_surf((j, i))
        
        # redo _find_surf_set() globally
        # cannot find a way to update it locally
        self._find_surf_set()
                    

    def plot(self, figsize=(8, 8), dpi=600, fname='Mesh.png'):
        """Plot mesh and surface."""
        fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi,
                                 constrained_layout=True)
        ax = axes[0]
        ax.scatter(self.x, self.z, c=self.mat, s=1, cmap=colMap, vmin=0.2)
        ax = axes[1]
        ax.scatter(self.x, self.z, c=self.surf, s=1)
        fig.savefig(fname, dpi=dpi)
        
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

    def change_mat(self, idx, mat_name='Vac'):
        """
        Update materials due to etching.

        idx: a.u., (i, j) tuple, the index of materails to be changed
        mat_name: str, material name
        """
        self.mat[idx] = mat_name

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
            sub_surf = np.concatenate((self.surf[bottom:top, left:], 
                           self.surf[bottom:top, 0:right]), axis=1)
            sub_x = np.concatenate((self.x[bottom:top, left:] - self.width, 
                                    self.x[bottom:top, 0:right]), 
                                   axis=1)
            sub_z = np.concatenate((self.z[bottom:top, left:], 
                                    self.z[bottom:top, 0:right]), axis=1)
        elif right > self.nx-1:
            right = right % self.nx
            sub_surf = np.concatenate((self.surf[bottom:top, left:], 
                           self.surf[bottom:top, 0:right]), axis=1)
            sub_x = np.concatenate((self.x[bottom:top, left:], 
                                    self.x[bottom:top, 0:right] + self.width), 
                                   axis=1)
            sub_z = np.concatenate((self.z[bottom:top, left:], 
                                    self.z[bottom:top, 0:right]), axis=1)
        else:
            sub_surf = self.surf[bottom:top, left:right]
            sub_x = self.x[bottom:top, left:right]
            sub_z = self.z[bottom:top, left:right]
        # print(sub_surf)

        if imode == "Fit Plane":
            # sub_surf consists of 1(surf) and -1(surf_vac)
            # surf_vac is not used when imode == 1, zero out -1
            temp_sub_surf = np.where(sub_surf == -1, 0, sub_surf)
            def cost_func_surf_norm(theta):
                """Construct the cost func for surface fitting."""
                A, B = -np.sin(theta), np.cos(theta)
                C = A*self.x[idx] + B*self.z[idx]
                Q = A*sub_x + B*sub_z - C
                Q = np.multiply(Q, temp_sub_surf)
                Q = np.power(Q, 2)
                Qsum = -abs(Q.sum())
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
            # sub_surf consists of 1(surf) and -1(surf_vac)
            # surf_vac is not used when imode == 2, zero out 1
            temp_sub_surf = np.where(sub_surf == 1, 0, sub_surf)
            # print(temp_sub_surf, sub_x)
            temp_vecx = np.multiply(self.x[idx] - sub_x, temp_sub_surf)
            temp_vecz = np.multiply(self.z[idx] - sub_z, temp_sub_surf)
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
                theta = np.arccos(surf_norm[0])
            else:
                theta = np.random.uniform(-np.pi, np.pi)
                surf_norm = np.array([sin(theta), -cos(theta)])
        
        return surf_norm, theta

    def plot_surf(self, figsize=(8, 8), dpi=600, fname='demo_surf.png',
                  surf_norm_range=2, surf_norm_mode='Fit Plane'):
        """Plot surf norm for all surf nodes."""
        fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi,
                                 constrained_layout=True)
        ax = axes[0]
        ax.scatter(self.x, self.z, c=self.mat, s=1, cmap=colMap, vmin=0.2)
        ax = axes[1]
        ax.scatter(self.x, self.z, c=self.surf, s=1)
        
        surf_list = np.transpose(np.nonzero(self.surf == 1))
        
        for surf_idx in surf_list:
            surf_idx = tuple(surf_idx)
            svec, sth = self.calc_surf_norm(surf_idx, radius=surf_norm_range, 
                                    imode=surf_norm_mode)        
            ax.quiver(self.x[surf_idx], self.z[surf_idx],
                          svec[0], svec[1], width=0.001)
           
        fig.savefig(fname, dpi=dpi)

    def find_float_cell(self, imode='Direct', idiag=0):
        """
        Search for the floating cells.
        
        Scan each material cell 
        imode: str, var, ['Direct', 'Diagonal']
        imode = 'Direct': check its 4 neighbours, top, bottom, left and right.
        imode = 'Diagonal': check additional 4 diagonal neighbours, top-left,
                            top-right, bottom-left, bottom-right.
        If they are all empty, the cell is identified as a floating cell, 
        which will be dropped to bottom.
        """
        float_cell = list()
        # if idiag == 1, record the floating cells
        if idiag:
            self.float_cell = np.zeros_like(self.mat).astype(int)
        for j in range(1, self.nz-1):
            for i in range(self.nx):
                # if mat[i,j] is not 0
                if self.mat[j, i]:
                    # check its 4 neighbours
                    temp = self.mat[j-1, i]
                    temp += self.mat[j+1, i]
                    temp += self.mat[j, i-1]
                    temp += self.mat[j, (i+1) % self.nx]
                    if not temp:
                        float_cell.append((j, i))
                        if idiag:
                            self.float_cell[j, i] = 1
        # plot the floating cells    
        if idiag:
            fig = plt.figure(figsize=(4, 8), dpi=300)
            plt.scatter(self.x, self.z, c=self.float_cell, s=1,
                        cmap=colMap, vmin=0.2)
            fig.savefig('float_cells.png', dpi=300)
            plt.show()

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
    
    def drop_floating_cell(self, imode='Remove'):
        """
        Drop the floating cells/clusters.
        
        imode: str, determine how to drop to the floating cells/clusters.
            imode = 'Remove', simply remove the floating cells/clusters.
            imode = 'Drop', drop the floating cells/clusthers downwards.
        """
        if imode == 'Remove':
            _idx_arr = np.where(self.surf == 1)
            pass

def rect_conv(coord, res_x, res_z):
    """Convert rectangular coordinates to index."""
    ncoord = [np.ceil(coord[0]/res_x), np.ceil(coord[1]/res_z),
              np.ceil(coord[2]/res_x), np.ceil(coord[3]/res_z)]
    ncoord = [int(temp) for temp in ncoord]
    return [ncoord[0], ncoord[0]+ncoord[2], ncoord[1], ncoord[1]+ncoord[3]]


if __name__ == '__main__':
    from FeatMod2d_ops import width, height, res_x, res_z
    from FeatMod2d_mat import mat0, mat1, mat2, Si2d
    mesh = MESHGRID(width, height, res_x, res_z)
    mesh.add_mat(Si2d)
    mesh.plot()
    
    colMap = copy.copy(cm.get_cmap("Accent"))
    colMap.set_under(color='white')

    rec_surf = []
    for temp_idx in mesh.surf_set:
        # temp_idx = (224, 13)
        temp_svec, temp_stheta = mesh.calc_surf_norm(temp_idx, 
                                                     radius=4,
                                                     imode='Sum Vector')
        rec_surf.append([temp_idx, temp_svec])

    def plot_surf_norm(ax, posn, svec):
        ax.quiver(posn[0], posn[1],
                  svec[0], svec[1])

    fig, axes = plt.subplots(1, 2, figsize=(8, 8),
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
