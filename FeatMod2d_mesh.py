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
        self.nx = math.ceil(width/res_x) # num of cells in x 
        self.nz = math.ceil(height/res_z) # num of cells iin z
        self.mat = np.zeros((self.nx, self.nz)) # mesh materials
        self.surf = np.array([]) # surface node
    
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
                self.mat[ncoord[0]:ncoord[1], 
                         ncoord[2]:ncoord[3]] = mater[1]
    
    def find_surface(self):
        surf = []
        # search surface within materials
        for i in range(self.nx):
            for j in range(1, self.nz-1):
                # if mat[i,j] is not 0
                if self.mat[i,j]:
                    # temp = multiply all neighbours
                    # temp = 0 means one of neighbours = 0
                    temp = self.mat[(i-1) % self.nx, j] \
                          *self.mat[(i+1) % self.nx, j] \
                          *self.mat[i, (j-1)% self.nz] \
                          *self.mat[i, (j+1) % self.nz]
                    if not temp:
                        surf.append((i,j))
        
        # search for left and right boundaries
        for i in [0, self.nx-1]:
            for j in range(1, self.nz-1):
                if not self.mat[i,j]:
                    surf.append((i,j))
        
        self.surf = np.array(surf).T
    
    def plot(self):        
        colMap = cm.Accent
        colMap.set_under(color='white')
        
        x = np.linspace(0.0, self.width, self.nx) 
        z = np.linspace(0.0, self.height, self.nz) 
        X, Z = np.meshgrid(x, z) 

        fig, axes = plt.subplots(1,2, figsize=(4,8),
                                   constrained_layout=True)
        axes[0].contourf(X, Z, self.mat.T, 
            cmap = colMap, vmin = 0.2, extend='both')
        axes[1].plot(self.surf[0, :]*self.res_x, self.surf[1, :]*self.res_z, 
            'o', )
        axes[1].set_xlim(0, self.width)
        axes[1].set_ylim(0, self.height)
        plt.show(fig)


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

