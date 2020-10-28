"""
Materials
"""

mat0 = []
m = [('SiO2', 1), 'rect', (0.0, 0.0, 100.0, 50.0)]
mat0.append(m)
m = [('Si', 2),   'rect', (0.0, 50.0, 100.0, 300.0)]
mat0.append(m)
m = [('PR', 3),   'rect', (0.0, 350.0, 30.0, 100.0)]
mat0.append(m)
m = [('PR', 3),   'rect', (70.0, 350.0, 30.0, 100.0)]
mat0.append(m)
m = [('Vac', 0),   'circ', (50.0, 350.0, 30.0)]
mat0.append(m)


mat1 = []
m = [('SiO2', 1), 'rect', (0.0, 0.0, 100.0, 50.0)]
mat1.append(m)
m = [('Si', 2),   'rect', (0.0, 50.0, 100.0, 300.0)]
mat1.append(m)
m = [('PR', 3),   'rect', (0.0, 350.0, 30.0, 50.0)]
mat1.append(m)
m = [('PR', 3),   'rect', (70.0, 350.0, 30.0, 50.0)]
mat1.append(m)
m = [('Vac', 0),   'circ', (50.0, 350.0, 30.0)]
mat1.append(m)

mat2 = []
m = [('SiO2', 1), 'rect', (0.0, 0.0, 100.0, 50.0)]
mat2.append(m)
m = [('Si', 2),   'rect', (0.0, 50.0, 100.0, 300.0)]
mat2.append(m)
m = [('PR', 3),   'rect', (0.0, 350.0, 30.0, 50.0)]
mat2.append(m)
m = [('PR', 3),   'rect', (70.0, 350.0, 30.0, 50.0)]
mat2.append(m)
m = [('Vac', 0),   'circ', (50.0, 350.0, 30.0)]
mat2.append(m)
m = [('Vac', 0),   'trgl', (30.0, 350.0, 30.0, 450.0, 0.0, 450.0)]
mat2.append(m)

# Si2d is designed to verify the code with MCFPM
Si2d = []
m = [('SiO2', 1), 'rect', (0.0, 0.0, 200.0, 25.0)]
Si2d.append(m)
m = [('Si', 2),   'rect', (0.0, 25.0, 200.0, 200.0)]
Si2d.append(m)
m = [('PR', 3),   'rect', (0.0, 225.0, 80.0, 25.0)]
Si2d.append(m)
m = [('PR', 3),   'rect', (120.0, 225.0, 80.0, 25.0)]
Si2d.append(m)
m = [('Vac', 0),   'trgl', (80.0, 225.0, 80.0, 250.0, 70.0, 250.0)]
Si2d.append(m)
m = [('Vac', 0),   'trgl', (120.0, 225.0, 120.0, 250.0, 130.0, 250.0)]
Si2d.append(m)
