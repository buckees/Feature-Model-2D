"""
Store species list.

Any species/particles used in the feature model should be selected from here.
"""

# name, ptype, mass, charge
Eon = {'name':'e', 'ptype':'Eon', 'mass':1.0, 'charge':-1}
H = {'name':'H', 'ptype':'Neut', 'mass':1.0, 'charge':0}
Hp = {'name':'H+', 'ptype':'Ion', 'mass':1.0, 'charge':1}
Ar = {'name':'Ar', 'ptype':'Neut', 'mass':32.0, 'charge':0}
Arp = {'name':'Ar+', 'ptype':'Ion', 'mass':32.0, 'charge':1}
He = {'name':'He', 'ptype':'Neut', 'mass':4.0, 'charge':0}
Hep = {'name':'He+', 'ptype':'Ion', 'mass':4.0, 'charge':1}
Cl2 = {'name':'Cl2', 'ptype':'Neut', 'mass':70.0, 'charge':0}
Cl2p = {'name':'Cl2+', 'ptype':'Ion', 'mass':70.0, 'charge':1}
Cl2n = {'name':'Cl2-', 'ptype':'Ion', 'mass':70.0, 'charge':-1}
Cl = {'name':'Cl', 'ptype':'Neut', 'mass':35.0, 'charge':0}
Clp = {'name':'Cl+', 'ptype':'Ion', 'mass':35.0, 'charge':1}
Cln = {'name':'Cl-', 'ptype':'Ion', 'mass':35.0, 'charge':-1}
O2 = {'name':'O2', 'ptype':'Neut', 'mass':32.0, 'charge':0}
O2p = {'name':'O2+', 'ptype':'Ion', 'mass':32.0, 'charge':1}
O2n = {'name':'O2-', 'ptype':'Ion', 'mass':32.0, 'charge':-1}
O = {'name':'O', 'ptype':'Neut', 'mass':16.0, 'charge':0}
Op = {'name':'O+', 'ptype':'Ion', 'mass':16.0, 'charge':1}
On = {'name':'O-', 'ptype':'Ion', 'mass':16.0, 'charge':-1}
F = {'name':'F', 'ptype':'Neut', 'mass':19.0, 'charge':0}
Si = {'name':'Si', 'ptype':'Neut', 'mass':28.0, 'charge':0}



sp_full_list = {
    'Eon':Eon,
    'H':H,
    'H+':Hp,
    'He':He,
    'He+':Hep,
    'Ar': Ar,
    'Ar+': Arp,
    'Cl2': Cl2,
    'Cl2+': Cl2p,
    'Cl2-': Cl2n,
    'Cl': Cl,
    'Cl+': Clp,
    'Cl-': Cln,
    'O2': O2,
    'O2+': O2p,
    'O2-': O2n,
    'O': O,
    'O+': Op,
    'O-': On,
    'F':F,
    'Si': Si
    }
