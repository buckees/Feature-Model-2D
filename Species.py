# -*- coding: utf-8 -*-
# Define the class of SP

class PARTICLE(object):
    """Stores all species"""
    def __init__(self, name, ptype, mass, charge):
        self.name = name # str
        self.ptype = ptype # str, 'E','Ion','Neut' or 'Bkg'
        self.mass = mass # unit in AMU
        self.charge = charge # unit in Unit Charge of Electron
    
    def __str__(self):
    return """
           This mesh with domain of %.1f nm (width) x %.1f nm (height)
           and resolution in (x, z) = (%.2f nm, %.2f nm)
           and number of cells in (x, z) = (%d, %d)
           """ \
            % (self.width, self.height, 
               self.res_x, self.res_z, 
               self.nx, self.nz)


# proton-electron mass ratio = 1836.15
Eon = SP('E',   'E',    5.45e-4, -1.0, 5.0) 
Arp = SP('Ar+', 'Ion',  32.0,     1.0, 0.2)
Ar  = SP('Ar',  'Bkg',  32.0,     0.0, 0.025)
Hp  = SP('H+',  'Ion',   1.0,     1.0, 0.2)
H   = SP('H',   'Bkg',   1.0,     0.0, 0.025)
N2  = SP('N2',  'Bkg',  28.0,     0.0, 0.025)