"""Feature Model 2D. Reflection."""

import numpy as np
from scipy.stats import cosine


class REFLECT(object):
    """Reflection Probability."""

    def __init__(self, ptcl='Ar+', mat='PR', prob=1.0):
        self.prob = prob  # reflection probability
        self.ptcl = ptcl  # particle species
        self.mat = mat  # materiasl which particle hits
        self.theta = 0.0  # reflection angle wrt surface normal

    def calc_prob(self):
        """Calc probability."""
        if self.mat == 'PR':
            self.prob = 1.0
        elif self.mat == 'Si':
            self.prob = 0.7
        return self.prob

    def rotate_random(self, uvec):
        """Rotate the direction randomly."""
        theta = np.random.uniform(-np.pi/4.0, + np.pi/4.0)
        theta = theta + np.pi
        x1, z1 = uvec
        x2 = np.cos(theta)*x1 - np.sin(theta)*z1
        z2 = np.sin(theta)*x1 + np.cos(theta)*z1
        return np.array([x2, z2])

    def diff_rflct(self, theta):
        """Calc the diffusive refelection."""
        self.theta  = cosine.rvs(size=1)
        self.theta += theta
        uvec = (np.cos(theta), np.sin(theta))
