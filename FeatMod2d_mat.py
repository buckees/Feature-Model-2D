"""
Materials
"""

mat0 = []
m = [('SiO2_', 1), 'rect', (0.0, 0.0, 100.0, 50.0)]
mat0.append(m)
m = [('Si_', 2),   'rect', (0.0, 50.0, 100.0, 300.0)]
mat0.append(m)
m = [('PR_', 3),   'rect', (0.0, 350.0, 30.0, 100.0)]
mat0.append(m)
m = [('PR_', 3),   'rect', (70.0, 350.0, 30.0, 100.0)]
mat0.append(m)
m = [('Vac_', 0),   'circ', (50.0, 350.0, 30.0)]
mat0.append(m)


mat1 = []
m = [('SiO2_', 1), 'rect', (0.0, 0.0, 100.0, 50.0)]
mat1.append(m)
m = [('Si_', 2),   'rect', (0.0, 50.0, 100.0, 300.0)]
mat1.append(m)
m = [('PR_', 3),   'rect', (0.0, 350.0, 30.0, 50.0)]
mat1.append(m)
m = [('PR_', 3),   'rect', (70.0, 350.0, 30.0, 50.0)]
mat1.append(m)
m = [('Vac_', 0),   'circ', (50.0, 350.0, 30.0)]
mat1.append(m)

mat2 = []
m = [('SiO2_', 1), 'rect', (0.0, 0.0, 100.0, 50.0)]
mat2.append(m)
m = [('Si_', 2),   'rect', (0.0, 50.0, 100.0, 300.0)]
mat2.append(m)
m = [('PR_', 3),   'rect', (0.0, 350.0, 30.0, 50.0)]
mat2.append(m)
m = [('PR_', 3),   'rect', (70.0, 350.0, 30.0, 50.0)]
mat2.append(m)
m = [('Vac_', 0),   'circ', (50.0, 350.0, 30.0)]
mat2.append(m)
m = [('Vac_', 0),   'trgl', (30.0, 350.0, 30.0, 450.0, 0.0, 450.0)]
mat2.append(m)

# Si2d is designed to verify the code with MCFPM
Si2d = []
m = [('SiO2_', 1), 'rect', (0.0, 0.0, 100.0, 25.0)]
Si2d.append(m)
m = [('Si_', 2),   'rect', (0.0, 25.0, 100.0, 200.0)]
Si2d.append(m)
m = [('PR_', 3),   'rect', (0.0, 225.0, 30.0, 25.0)]
Si2d.append(m)
m = [('PR_', 3),   'rect', (70.0, 225.0, 30.0, 25.0)]
Si2d.append(m)
m = [('Vac_', 0),   'trgl', (30.0, 225.0, 30.0, 250.0, 25.0, 250.0)]
Si2d.append(m)
m = [('Vac_', 0),   'trgl', (70.0, 225.0, 70.0, 250.0, 75.0, 250.0)]
Si2d.append(m)

# Si2d_trench is designed to verify the bottom reflection
Si2d_trench = Si2d.copy()
m = [('Vac_', 0),   'rect', (30.0, 25.0, 40.0, 200.0)]
Si2d_trench.append(m)

# Si2d_trench is designed to verify the bottom reflection
Si2d_trench_v02 = Si2d.copy()
m = [('Vac_', 0),   'rect', (30.0, 100.0, 40.0, 200.0)]
Si2d_trench_v02.append(m)
