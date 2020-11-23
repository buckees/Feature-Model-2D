"""Feature Model 2D. Partile file."""

import numpy as np
from math import cos, sin, sqrt, acos
import matplotlib.pyplot as plt
from scipy.stats import cosine


class PARTICLE(object):
    """Create particle object."""

    def __init__(self, name, ptype, mass, charge, dead=0):
        self.name = name  # str
        self.ptype = ptype  # str, 'E','Ion','Neut' or 'Bkg'
        self.mass = mass  # unit in AMU
        self.charge = charge  # unit in Unit Charge of Electron
        self.posn = np.zeros(2)
        self.enrg = 0.025  # unit in eV, initial as room temperature
        self.uvec = np.zeros(2)
        self.accl = np.zeros(2)
        self.dead = dead  # indicator for ptcl alive or dead

    def __str__(self):
        """Print out the particle informaiton."""
        return """
               The launched particle is '%s' with mass of %.1f AMU
               it belongs to '%s', with charge of %+d
               """ \
               % (self.name, self.mass,
                  self.ptype, self.charge)

    def init_posn(self, width, height):
        """Initialize the position at the top boundary."""
        init_posx = np.random.uniform(0.0, 1.0)*width
        self.posn = np.array([init_posx, height])

    def init_uvec(self, idstrb=None):
        """
        Initialize the velocity direction.

        in (x, -z) (half-down quadrant).
        """
        itype = idstrb[0]
        if itype == 'Uniform2D':
            left, right = idstrb[1]/180.0*np.pi, idstrb[2]/180.0*np.pi
            theta = np.random.uniform(left, right)
        elif itype == 'Uniform3D':
            left, right = idstrb[1]/180.0*np.pi, idstrb[2]/180.0*np.pi
            temp_th = np.random.uniform(left, right)
            temp_phi = np.random.uniform(-np.pi, np.pi)
            temp_vz = cos(temp_th)
            temp_vx = sin(temp_th)*cos(temp_phi)
            # temp_vy = sin(temp_th)*sin(temp_phi)
            temp_v = sqrt(temp_vz**2 + temp_vx**2)
            temp_vz = temp_vz/temp_v
            theta = acos(temp_vz)*np.sign(sin(temp_phi))
        elif itype == 'Normal':
            mu, sigma = 0.0, 0.1  # default mean and standard deviation
            mu, sigma = idstrb[1], idstrb[2]
            theta = np.random.normal(mu, sigma)
        elif itype == 'Cosine':
            scale = idstrb[1]/180.0*np.pi
            theta = cosine.rvs(scale=scale, size=1)
        elif itype == 'Mono':
            theta = idstrb[1]/180.0*np.pi

        self.uvec = np.array([sin(theta), -cos(theta)])

    def init_enrg(self, idstrb='Uniform',
                  enrg_min=1e-2, enrg_max=1e4):
        """Initialize the particle energy."""
        if idstrb == 'Uniform':
            enrg = np.random.uniform(enrg_min, enrg_max)

        elif idstrb == 'Normal':
            pass

        elif idstrb == 'Cosine':
            mu, sigma = 0, 0.1  # mean and standard deviation
            pass

        self.enrg = enrg

    def init_plot(self):
        """Plot the initialized positions and velocities."""
        fig, axes = plt.subplots(1, 2, figsize=(10, 4),
                                 constrained_layout=True)
        axes[0].plot(self.uvec[0], self.uvec[1], 'o')
        axes[0].set_xlim(-1.0, 1.0)
        axes[0].set_ylim(-1.0, 0.0)
        axes[1].quiver(self.posn[0], self.posn[1],
                       self.uvec[0], self.uvec[1])
        plt.show()

    def move_ptcl(self, delta_L):
        """Move each partile in a length of delta_L along its v-vector."""
        self.posn += self.uvec*delta_L

    def bdry_check(self, width, height, imode='lost'):
        """
        Check the b.c. for moving ptcl.

        make the ptcl dead if it gets beyond the top bdry
        three modes for vertical bdry are available:
            lost, periodic and reflective
        left bdry is alway at 0.0
        """
        if self.posn[1] >= height:
            self.dead = 1
        elif imode == 'lost':
            if not (0.0 < self.posn[0] < width):
                self.dead = 1
        elif imode == 'periodic':
            self.posn[0] = self.posn[0] % width
        elif imode == 'reflective':
            pass


if __name__ == '__main__':
    from Species import Arp
    fig = plt.figure()
    for i in range(100):
        Arp.init_uvec(['Uniform3D', -45.0, 45.0])
        vec = Arp.uvec
        plt.scatter(vec[0], vec[1])
    plt.show()
        
    
