"""
Store species list.

Any species/particles used in the feature model should be selected from here.
"""

from FeatMod2d_ptcl import PARTICLE

# Ar = PARTICLE('Ar', 'Neut',  32.0,     0)
# Arp = PARTICLE('Ar+', 'Ion',  32.0,     1)
# Cl2 = PARTICLE('Cl2', 'Neut',  70.0,     0)
# Cl = PARTICLE('Cl', 'Neut',  35.0,     0)
# Clp = PARTICLE('Cl+', 'Ion',  35.0,     1)

# name, ptype, mass, charge
Ar = {'name':'Ar', 'ptype':'Neut', 'mass':32.0, 'charge':0}
Arp = {'name':'Ar+', 'ptype':'Ion', 'mass':32.0, 'charge':1}
Cl2 = {'name':'Cl2', 'ptype':'Neut', 'mass':70.0, 'charge':0}
Cl = {'name':'Cl', 'ptype':'Neut', 'mass':35.0, 'charge':0}
Clp = {'name':'Cl+', 'ptype':'Ion', 'mass':35.0, 'charge':1}

sp_full_list = {'Ar': Ar,
                'Ar+': Arp,
                'Cl2': Cl2,
                'Cl': Cl,
                'Cl+': Clp}

Si = PARTICLE('Si', 'Neut',  28.0,     0)