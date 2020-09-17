# -*- coding: utf-8 -*-
"""
Feature Model 2D
Particle Initiation
"""
import numpy as np
import math

from FeatMod2d_ops import *

def random_2d_vector(num_ptcl):
    """
    Generates a random 2D unit vector (direction) 
    with a uniform spherical distribution
    """
    theta = np.random.uniform(-np.pi/2.0, np.pi/2.0, size=num_ptcl)
    x = np.sin(theta)
    z = -np.cos(theta)
    return (x, z)

def random_3d_vector():
    """
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return (x,y,z)