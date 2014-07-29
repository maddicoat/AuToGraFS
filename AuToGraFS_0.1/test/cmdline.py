import numpy as np
from ase.tasks.main import run
atoms, task = run('H2 --bond-length=0.78 -R 0.01 -F 5,2 --atomize')
atoms, task = run('H2 H -s')
results = np.array(task.results['H2'])
assert abs(results -
           [1.071, 0.000, 0.779, 862.780, 5.349, 5.349]).max() < 0.001

atoms, task = run('bulk Cu -F 5,2')
atoms, task = run('bulk Cu -s')
results = np.array(task.results['Cu'])
assert abs(results -
           [-0.0057, 0.0014, 11.5654, 134.4389]).max() < 0.001
