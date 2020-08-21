# -*- coding: utf-8 -*-
"""
Feature Model 2D
Main program
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from FeatMod2d_ops import width, height, res_x, res_z, num_ptcl
from FeatMod2d_mesh import MESHGRID
from FeatMod2d_init import random_2d_vector
from FeatMod2d_move import find_intersect_node
from FeatMod2d_plot import plot_mesh, plot_itsct

# create mesh
mesh = MESHGRID(width, height, res_x, res_z)
mesh.mat_input()
mesh.find_surface()
mesh.plot()

