# -*- coding: utf-8 -*-
"""
Feature Model 2D
Reflection
"""

class REFLECT(object):
    """Reflection Probability"""
    def __init__(self, ptcl='Ar+', mat='PR', prob=1.0):
        self.prob = prob # reflection probability
        self.ptcl = ptcl # particle species
        self.mat = mat # materiasl which particle hits
        
    def calc_prob(self):
        if self.mat == 'PR':
            self.prob = 1.0
        elif self.mat == 'Si':
            self.prob = 0.1
        return self.prob
        
    
