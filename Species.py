"""
Store species list.

Any species/particles used in the feature model should be selected from here.
"""

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