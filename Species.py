"""
Store species list.

Any species/particles used in the feature model should be selected from here.
"""

from FeatMod2d_ptcl import PARTICLE

Ar = PARTICLE('Ar', 'Neut',  32.0,     0)
Arp = PARTICLE('Ar+', 'Ion',  32.0,     1)
Cl2 = PARTICLE('Cl2', 'Neut',  70.0,     0)
Cl = PARTICLE('Cl', 'Neut',  35.0,     0)
Clp = PARTICLE('Cl+', 'Ion',  35.0,     1)

Si = PARTICLE('Si', 'Neut',  28.0,     0)