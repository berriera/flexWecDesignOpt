import numpy as np
import pygmsh
import os
from analysis import boundary_condition_frequency_solver
from analysis import run_wamit
from file_mgmt import create_case_directory
from substitution import create_case_files
from mesh import create_mesh_file
from mesh import submerged_mesh

# Samples for Monte Carlo random sampling
N = 30

# File management
device_name = 'FlexibleBarge'
run_wamit_command = 'C:\WAMITv7\wamit'
common_file_directory = os.path.abspath(os.path.join('examples', 'barge', 'pygmsh_meshing'))
output_directory = os.path.abspath(os.path.join('examples', 'output'))

# Optional meshing arguments and WAMIT command
gmsh_exe_location = 'C:/Users/13365/Documents/gmsh-4.5.6-Windows64/gmsh-4.5.6-Windows64/gmsh'
mesh_refinement_factor = 0.40
run_wamit_command = 'C:\WAMITv7\wamit'

# File path to absolute location
common_file_directory = os.path.abspath(common_file_directory)
output_directory = os.path.abspath(output_directory)


# Roots of beam equation with free-free boundary conditions
def characteristic_equation(x):
    import math
    return math.cos(x) * math.cosh(x) - 1


# Finds 8 roots of beam equation
eigenvalues = boundary_condition_frequency_solver(characteristic_equation, 8)


class FlexibleBarge(object):

    def __init__(self, design_vars):
        self.device_name = 'barge'
        self.L = design_vars[0]
        self.w = design_vars[1]
        self.h = design_vars[2]

        # Barge material properties
        self.E = 30.720e6  # Modulus of elasticity, Pa
        self.rho = 500  # density, kg / (m ** 3)
        self.nu = 0.3

        # Barge geometry properties
        self.inertia_area = (1 / 12) * self.w * self.h ** 3
        self.cross_section_area = self.h * self.w

        # Barge inertial properties
        self.mass = self.rho * self.L * self.w * self.h
        self.Cg = 0
        self.kx = ((1 / 12) * (self.L ** 2 + self.h ** 2)) ** (1 / 2)
        self.ky = ((1 / 12) * (self.w ** 2 + self.h ** 2)) ** (1 / 2)
        self.kz = ((1 / 12) * (self.L ** 2 + self.w ** 2)) ** (1 / 2)
        self.Ix = self.mass * self.kx ** 2
        self.Iy = self.mass * self.ky ** 2
        self.Iz = self.mass * self.kz ** 2

        # Barge structural deformation mass M and stiffness C matrices
        self.M = self.mass * np.asarray([0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25])
        self.eigenvalues = eigenvalues
        self.kappa = self.eigenvalues / self.L
        self.omega = (((self.E * self.inertia_area) / (self.rho * self.cross_section_area * self.L ** 4)) *
                      self.eigenvalues ** 4) ** (1 / 2)
        self.C = (self.omega ** 2) * self.M
        # TODO: figure out where last value M(5,5) in barge.frc comes from;
        # Update: I think this value is just wrong and should be M(5,5) = Iy = 2.167E+09

    def substitutions(self):
        # Dictionary of values
        return {'L': self.L, 'w': self.w, 'h': self.h, 'Cg': self.Cg,
                'mass': self.mass,
                'M11': self.M[0], 'M22': self.M[1], 'M33': self.M[2], 'M44': self.M[3],
                'M55': self.M[4], 'M66': self.M[5], 'M77': self.M[6], 'M88': self.M[7],
                'C11': self.C[0], 'C22': self.C[1], 'C33': self.C[2], 'C44': self.C[3],
                'C55': self.C[4], 'C66': self.C[5], 'C77': self.C[6], 'C88': self.C[7]
                }

    def geometry(self):
        # Mesh properties
        geometry = pygmsh.opencascade.Geometry()
        geometry.add_box(x0=[-1 / 2 * self.L, -1 / 2 * self.w, -1 / 2 * self.h],
                         extents=[self.L, self.w, self.h])
        return geometry


# Monte Carlo sampling with seed for reproducibility
np.random.seed(42)
lower_bounds = np.asarray([10, 3, 3])
upper_bounds = np.asarray([200, 30, 30])

for i in range(N):
    print('Case: ', str(i + 1))
    design_variables = lower_bounds + np.random.random(size=3, ) * upper_bounds
    barge = FlexibleBarge(design_variables)
    barge_substitutions = barge.substitutions()
    barge_geometry = barge.geometry()
    create_case_directory(output_directory, i + 1)
    create_case_files(common_file_directory, barge_substitutions)
    # create_mesh_file(barge_geometry, device_name, gmsh_exe_location, mesh_refinement_factor)
    # submerged_mesh(device_name)
    # run_wamit(run_wamit_command)
print('Done.')
