"""
Feature Model 2D
Mesh file
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class MESHGRID(object):
    """Mesh object"""
    def __init__(self, width, height, res_x, res_z):
        self.width = width # domain width in x direction
        self.height = height # domain height in z direction
        self.res_x = res_x # resolution in x direction
        self.res_z = res_z # resolution in z direction
        self.res = np.array([self.res_x, self.res_z]) # array of resolution
#        <--------------------- width --------------------->
#        0.0 ---- 1.0 ---- 2.0 ---- ... ---- 99.0 ---- 100.0
#        --- cell --- cell --- cell ... cell ---- cell -----
#        -- center - center - center . center -- center ----
#        --- 0.5 ---- 1.5 ---- 2.5  ... 98.5 ---- 99.5 -----
#        --(0, nj)--(1, nj)--(2, nj)...(98, nj)--(99, nj)---
        self.nx = math.ceil(width/res_x) # num of cells in x 
        self.nz = math.ceil(height/res_z) # num of cels iin z
        # init x and z coordinates
        tempx = np.linspace(0.0, self.width, self.nx) 
        tempz = np.linspace(0.0, self.height, self.nz) 
        self.x, self.z = np.meshgrid(tempx, tempz)
        # mesh materials is assigned to self.mat matrix
        # note that the shape of self.mat is (nz, nx)
        self.mat = np.zeros_like(self.x).astype(int)
        self.surf = np.array([]) # surface node
        self.mater = [] # materials name <--> materails No.
    
    def __str__(self):
        return """
               This mesh with domain of %.1f nm (width) x %.1f nm (height)
               and resolution in (x, z) = (%.2f nm, %.2f nm)
               and number of cells in (x, z) = (%d, %d)
               """ \
                % (self.width, self.height, 
                   self.res_x, self.res_z, 
                   self.nx, self.nz)
    
    # assign input materials to the mesh
    def mat_input(self):
        self.mater = ['Vac','SiO2', 'Si', 'PR']
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
    
    def find_surface(self):
        surf = []
        # search surface within materials
        for j in range(1, self.nz-1):
            for i in range(self.nx):
                # if mat[i,j] is not 0
                if self.mat[j,i]:
                    # temp = multiply all neighbours
                    # temp = 0 means one of neighbours = 0
                    temp = self.mat[j, (i-1) % self.nx] \
                          *self.mat[j, (i+1) % self.nx] \
                          *self.mat[(j-1) % self.nz, i] \
                          *self.mat[(j+1) % self.nz, i]
                    if not temp:
                        surf.append((j,i))
        
        # search for left and right boundaries
        for j in range(1, self.nz-1):
            for i in [0, self.nx-1]:            
                if not self.mat[j,i]:
                    surf.append((j,i))
        
        self.surf = np.array(surf).T
    
    def plot(self):        
        colMap = cm.Accent
        colMap.set_under(color='white')
        
        fig, axes = plt.subplots(1,2, figsize=(4,8),
                                   constrained_layout=True)
        axes[0].contourf(self.x, self.z, self.mat, 
            cmap = colMap, vmin = 0.2, extend='both')
        axes[1].plot(self.surf[0, :]*self.res_x, self.surf[1, :]*self.res_z, 
            'o', )
        axes[1].set_xlim(0, self.width)
        axes[1].set_ylim(0, self.height)
        plt.show(fig)

    def hit_check(self, posn):
        """
        all posn needs to be rounded to 0.5*nx (cell center)
        posn needs to be shifted by half res_x
        """
        idx = np.rint((posn-self.res_x*0.5)/self.res).astype(int)
        idx = tuple(idx)
        return self.mat[idx], idx

    def update_surf(self, idx, threshold):
        """
        update surface due to etching
        threshold: etching probability
        mat == 2 equivalent to mat == 'Si'
        """
        if self.mat[idx] == 2:
                rnd = np.random.uniform(0.0, 1.0)
                if rnd < threshold:
                    self.mat[idx] = 0
    
    def surf_norm(self, idx):
        vec_norm = (0.0, 1.0)
        return vec_norm


# convert rectangular coordinates to index
def rect_conv(coord, res_x, res_z):
    ncoord = [math.ceil(coord[0]/res_x), math.ceil(coord[1]/res_z), 
              math.ceil(coord[2]/res_x), math.ceil(coord[3]/res_z)]
    return [ncoord[0], ncoord[0]+ncoord[2], ncoord[1], ncoord[1]+ncoord[3]]
        
if __name__ == '__main__':
    from FeatMod2d_ops import width, height, res_x, res_z
    mesh = MESHGRID(width, height, res_x, res_z)
    mesh.mat_input()
    mesh.find_surface()
    mesh.plot()

