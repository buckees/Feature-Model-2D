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

    def rand_rflct(self):
        """
        Calc the uniform random refelection.
        
        stheta = surface normal theta
        rtheta = reflective theta wrt x=0+
        return reflective unit vector
        """
        rtheta = np.random.uniform(-np.pi/2.0, np.pi/2.0)
        rtheta += self.stheta
        return np.array([np.cos(rtheta), np.sin(rtheta)])
    
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

    def mix_rflct(self, ivec, ratio=1.0):
        """
        Calc the mixed reflection from specular and diffusive.
        
        ivec = incident vector
        spec_uvec = specular reflective vector
        diff_uvec = diffusive reflective vector
        mix_uvec = mixed reflective vector
        """
        spec_uvec = self.spec_rflct(ivec)
        diff_uvec = self.diff_rflct()
        mix_uvec = ratio/(ratio+1.0)*spec_uvec + 1.0/(ratio+1.0)*diff_uvec
        return mix_uvec/np.linalg.norm(mix_uvec)
        
        