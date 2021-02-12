import numpy as np
import os

from flexWecDesignOpt.device import Device
from examples.barge.flexible_barge_class import FlexibleBarge

# Monte Carlo sampling with seed for reproducibility
sample_count = 30
design_variable_count = 5
objectives_count = 1
lower_bounds = np.asarray([10, 3, 3, 20e6, 250])
upper_bounds = np.asarray([200, 30, 30, 40e6, 750])
np.random.seed(42)

# Set up data collection
design_variables = np.zeros(shape=(sample_count, design_variable_count))
f = np.zeros(shape=(sample_count, objectives_count))

# Set constant class variables
Device.name = 'barge'
Device.input_file_location = os.path.abspath(r'input_files')
Device.copy_file_list = ['fnames.wam', 'config.wam', 'barge.cfg']
Device.edit_file_list = ['barge.frc', 'barge_Length.dat']
Device.output_file_location = r'C:\Users\13365\Desktop\optimization_output'
Device.gmsh_location = r'C:\Users\13365\Documents\gmsh-4.5.6-Windows64\gmsh-4.5.6-Windows64\gmsh.exe'

for i in range(sample_count):
    print('Evaluating design {}:'.format(i + 1))
    design_variables[i, :] = lower_bounds \
                             + np.random.random(size=(design_variable_count,)) * (upper_bounds - lower_bounds)
    barge = FlexibleBarge(design_variables[i, :])

    device = Device(design_variables)
    device.mesh_geometry = barge.barge_geometry()
    device.create_case_directory()
    device.create_case_files(barge.barge_substitutions())
    device.mesh()
    # device.run_wamit()

    f[i:, 0] = barge.calculate_objective()  # TODO: finish objective function parameters
print('\n\nDone with all designs.')
f_min = np.min(f)
design_min = design_variables[np.where(f == f_min)[0][0], :]

print('\n\tObjective function minimum: {}'.format(np.min(f)))
print('Design variables corresponding to minimum: ')
print(design_min)

# TODO: add a .csv import example somewhere
