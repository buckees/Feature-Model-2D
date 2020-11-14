"""Feature Model 2D. Reflection."""

import numpy as np
from math import cos, sin
from scipy.stats import cosine


class REFLECT(object):
    """Reflection Probability."""

    def __init__(self, ptcl='Ar+', mat='PR', prob=1.0, 
                 svec=np.array([0.0, 1.0]), stheta=np.pi/2.0):
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
            self.prob = 0.5
        return self.prob

    def rflct(self, ivec, erg=1.0, imode={'spec':0.3, 'mix':0.4 , 'diff':0.3}):
        """
        Calc the reflection.
        
        ivec: a.u., (2, ) array, incident vector
        erg: eV, var, incident particle energy
        imode: a.u., dict, determine the prob of each rflct mode
        return a vec: a.u., (2, ) array, reflected vec
        """
        rnd = np.random.uniform(0.0, 1.0)
        if 0.0 <= rnd < imode['spec']:
            return self._spec_rflct(ivec)
        elif  imode['spec'] <= rnd < (imode['spec'] + imode['mix']):
            return self._mix_rflct(ivec, ratio=(0.7-rnd)/(rnd-0.3) )
        else:
            return self._diff_rflct()

    def _revs_rflct(self, ivec):
        """Reverse the incident vector."""
        return -ivec

    def _rand_rflct(self):
        """
        Calc the uniform random refelection.
        
        stheta = surface normal theta
        rtheta = reflective theta wrt x=0+
        return reflective unit vector
        """
        rtheta = np.random.uniform(-np.pi/2.0, np.pi/2.0)
        rtheta += self.stheta
        return np.array([np.cos(rtheta), np.sin(rtheta)])
    
    def _spec_rflct(self, ivec):
        """
        Calc the specular reflection.
        
        svec: surface normal vector.
        rvec: reflective vector 
        ivec: incident vector
        rvec = ivec - 2*(ivec dot svec)*svec
        """
        return ivec - 2.0*np.dot(ivec, self.svec)*self.svec
        
    def _diff_rflct(self):
        """
        Calc the diffusive refelection.
        
        stheta = surface normal theta
        rtheta = reflective theta wrt x=0+
        return reflective unit vector
        """
        rtheta, = cosine.rvs(size=1, scale=0.5)
        rtheta += self.stheta
        return np.array([cos(rtheta), sin(rtheta)])

    def _mix_rflct(self, ivec, ratio=1.0):
        """
        Calc the mixed reflection from specular and diffusive.
        
        ivec = incident vector
        ratio = ratio of specular to diffusive
        spec_uvec = specular reflective vector
        diff_uvec = diffusive reflective vector
        mix_uvec = mixed reflective vector
        """
        spec_uvec = self._spec_rflct(ivec)
        diff_uvec = self._diff_rflct()
        mix_uvec = ratio/(ratio+1.0)*spec_uvec + 1.0/(ratio+1.0)*diff_uvec
        return mix_uvec/np.linalg.norm(mix_uvec)
        
        