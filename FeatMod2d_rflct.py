"""Feature Model 2D. Reflection."""

import numpy as np
from math import pi, sin, cos


class REFLECT(object):
    """Reflection Probability."""

    def __init__(self, ptcl='Ar+', mat='PR', prob=1.0):
        self.prob = prob  # reflection probability
        self.ptcl = ptcl  # particle species
        self.mat = mat  # materiasl which particle hits

    def calc_prob(self):
        """Calc probability."""
        if self.mat == 'PR':
            self.prob = 1.0
        elif self.mat == 'Si':
            self.prob = 0.1
        return self.prob

    def rotate_random(self, uvec):
        """Rotate the direction randomly."""
        theta = np.random.uniform(-pi/4.0, + pi/4.0)
        theta = theta + pi
        x1, z1 = uvec
        x2 = cos(theta)*x1 - sin(theta)*z1
        z2 = sin(theta)*x1 + cos(theta)*z1
        return np.array([x2, z2])
