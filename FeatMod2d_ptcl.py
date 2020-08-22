"""
Feature Model 2D
Partile file
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class PARTICLE(object):
    """Stores all species"""
    def __init__(self, name, ptype, mass, charge):
        self.name = name # str
        self.ptype = ptype # str, 'E','Ion','Neut' or 'Bkg'
        self.mass = mass # unit in AMU
        self.charge = charge # unit in Unit Charge of Electron
        self.posn = np.zeros(2)
        self.vels = np.zeros(2)
        self.uvec = np.zeros(2)
        self.accl = np.zeros(2)
    
    def __str__(self):
        return """
               The launched particle is '%s' with mass of %.1f AMU
               it belongs to '%s', with charge of %+d
               """ \
               % (self.name, self.mass, 
                  self.ptype, self.charge)
    
    def init_posn(self, width, height):
        """
        Initialize the position at the top boundary
        """
        init_posx = np.random.uniform(0.0, 1.0)*width
        self.posn = np.array([init_posx, height])
    
    def init_vels(self):
        """
        Initialize the velocity in (x, -z) (half-down quadrant)
        """
        theta = np.random.uniform(-np.pi/2.0, np.pi/2.0)
        self.uvec = np.array([math.sin(theta), -math.cos(theta)])
    
    def init_plot(self):
        """
        Plot the initialized positions and velocities
        """
        fig, axes = plt.subplots(1,2, figsize=(10,4),
                                   constrained_layout=True)
        axes[0].plot(self.uvec[0], self.uvec[1], 'o') 
        axes[0].set_xlim(-1.0, 1.0)
        axes[0].set_ylim(-1.0, 0.0)
        axes[1].quiver(self.posn[0], self.posn[1],
                       self.uvec[0], self.uvec[1])
        plt.show(fig)

    def move_ptcl(self, delta_L):
        """
        Move each partile in a length of delta_L along its v-vector
        """
        self.posn += self.uvec*delta_L

if __name__ == '__main__':
    from FeatMod2d_ops import width, height, res_x, res_z
    Arp = PARTICLE('Ar+', 'Ion',  32.0,     1)
    print(Arp)
    Arp.init_posn(width, height)
    Arp.init_vels()
    Arp.init_plot()
    
    delta_L = min(res_x, res_z)
    for i in range(100):
        Arp.move_ptcl(delta_L)
    Arp.init_plot()
    
