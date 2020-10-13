"""Feature Model 2D. Reflection."""

import numpy as np
from scipy.stats import cosine


class REFLECT(object):
    """Reflection Probability."""

    def __init__(self, ptcl='Ar+', mat='PR', prob=1.0):
        self.prob = prob  # reflection probability
        self.ptcl = ptcl  # particle species
        self.mat = mat  # materiasl which particle hits
        self.svec = svec  # surface normal vector at hit
        self.stheta = stheta  # surface normal angle at hit

    def calc_prob(self):
        """Calc probability."""
        if self.mat == 'PR':
            self.prob = 1.0
        elif self.mat == 'Si':
            self.prob = 0.7
        return self.prob

    def rotate_random(uvec):
        """Rotate the direction randomly."""
        theta = np.random.uniform(-np.pi/4.0, + np.pi/4.0)
        theta = theta + np.pi
        x1, z1 = uvec
        x2 = np.cos(theta)*x1 - np.sin(theta)*z1
        z2 = np.sin(theta)*x1 + np.cos(theta)*z1
        return np.array([x2, z2])
    
    def spec_rflct(self, ivec):
        """
        Calc the specular reflection.
        
        svec: surface normal vector.
        rvec: reflective vector 
        ivec: incident vector
        rvec = ivec - 2*(ivec dot svec)*svec
        """
        return ivec - 2.0*np.dot(ivec, self.svec)*self.svec
        
    def diff_rflct(self):
        """
        Calc the diffusive refelection.
        
        stheta = surface normal theta
        rtheta = reflective theta wrt x=0+
        return reflective unit vector
        """
        rtheta = cosine.rvs(size=1)[0]
        rtheta += self.stheta
        return np.array([np.cos(rtheta), np.sin(rtheta)])
