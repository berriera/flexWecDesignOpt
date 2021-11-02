import numpy as np
import os
import pyswarms as ps

from flexWecDesignOpt.device import Device
from examples.barge.flexible_barge_class import FlexibleBarge

lower_bounds = np.asarray([10, 10, 10, 20e6, 250])
upper_bounds = np.asarray([20, 30, 30, 40e6, 750])
bounds = (lower_bounds, upper_bounds)
np.random.seed(42)

# Set constant class variables
Device.name = 'barge'
Device.input_file_location = os.path.abspath(r'input_files')
Device.copy_file_list = ['fnames.wam', 'config.wam', 'barge.cfg']
Device.edit_file_list = ['barge.frc', 'barge_Length.dat']
Device.output_file_location = r'C:\Users\13365\Desktop\optimization_output'
Device.gmsh_location = r'C:\Users\13365\Documents\gmsh-4.5.6-Windows64\gmsh-4.5.6-Windows64\gmsh.exe'


def evaluate_barge(design_variables):
    n_particles = design_variables.shape[0]
    f = []

    for k in range(n_particles):
        particle_variables = design_variables[k, :]
        evaluate_barge.counter += 1
        barge = FlexibleBarge(particle_variables)

        device = Device(particle_variables)
        # device.mesh_geometry = barge.barge_geometry()
        # device.create_case_directory()
        # device.create_case_files(barge.barge_substitutions())
        # device.mesh()
        # device.run_wamit()
        f.append(barge.calculate_objective())

    return np.array(f)  # TODO: finish objective function parameters


evaluate_barge.counter = 0
options = {'c1': 0.8, 'c2': 0.3, 'w': 0.5}
optimizer = ps.single.GlobalBestPSO(n_particles=60, dimensions=5, options=options, bounds=bounds)

best_cost, best_position = optimizer.optimize(evaluate_barge, iters=1000)

print('\nObjective function minimum: {}'.format(best_cost))
print('\nDesign variables corresponding to minimum: ')
print(best_position)

print("Iteration count: {}".format(evaluate_barge.counter))
