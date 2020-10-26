"""
Materials
"""

materials = []
m = [('SiO2', 1), 'rect', (0.0, 0.0, 100.0, 50.0)]
materials.append(m)
m = [('Si', 2),   'rect', (0.0, 50.0, 100.0, 300.0)]
materials.append(m)
m = [('PR', 3),   'rect', (0.0, 350.0, 30.0, 100.0)]
materials.append(m)
m = [('PR', 3),   'rect', (70.0, 350.0, 30.0, 100.0)]
materials.append(m)
m = [('Vac', 0),   'circ', (50.0, 350.0, 30.0)]
materials.append(m)