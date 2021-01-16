import numpy as np
from examples.barge.flexible_barge_class import FlexibleBarge
#from flexible_barge_example import FlexibleBarge
import analysis

# Monte Carlo sampling with seed for reproducibility
sample_count = 30
design_variable_count = 5
objectives_count = 1
lower_bounds = np.asarray([10, 3, 3, 20e6, 250])
upper_bounds = np.asarray([200, 30, 30, 40e6, 750])
np.random.seed(42)

# Set up
design_variables = np.zeros(shape=(sample_count, design_variable_count))
f = np.zeros(shape=(sample_count, objectives_count))

for i in range(sample_count):
    print('Evaluating design {} with variables'.format(i + 1))
    design_variables[i, :] = lower_bounds + np.random.random(size=(design_variable_count, )) * (upper_bounds - lower_bounds)
    print(design_variables[i, :])
    f[i, :] = analysis.evaluate_device(device_class=FlexibleBarge, design_variables=design_variables[i, :])
print('Done.')
# TODO: print min objective function with design variables

# Results
print('Objective function values: ')
print(f)
print('Design variables: ')
print(design_variables)

# TODO: add a .csv import example somewhere
