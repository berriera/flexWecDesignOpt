import pygmsh
import os
import numpy as np
from file_mgmt import create_case_directory
from substitution import create_case_files
from mesh import create_mesh_file
from mesh import submerged_mesh
from analysis import run_wamit

run_wamit_command = 'C:\WAMITv7\wamit'
gmsh_exe_location = 'C:/Users/13365/Documents/gmsh-4.5.6-Windows64/gmsh-4.5.6-Windows64/gmsh'
common_file_directory = os.path.abspath(os.path.join('examples', 'flexible_tube'))
output_directory = os.path.abspath(os.path.join('examples', 'output'))
device_name = 'tube'
mesh_refinement_factor = 0.50


class FlexibleTube(object):

    def __init__(self, design_vars):
        # Simulation parameters
        self.seafloor_depth = 5.0
        self.resonant_mode_count = 1
        self.degrees_of_freedom = 6 + self.resonant_mode_count

        # Design variable array
        self.device_name = 'tube'
        self.depth = design_vars[0]
        self.radius = design_vars[1]
        self.thickness = design_vars[2]
        self.length = design_vars[3]
        self.distensibility = design_vars[4]
        self.eta = design_vars[5]
        self.towhead_mass = design_vars[6]

        # Mass properties
        self.rho = 532.6462876469635  # Calculated from Babarit et al. 2017
        self.tube_volume = 2 * self.rho * np.pi * self.radius * self.thickness * self.length
        self.tube_mass = self.rho * self.tube_volume
        self.mass = self.tube_mass + 2 * self.towhead_mass
        self.Ix = self.tube_mass * (self.radius ** 2) + 2 * ((1 / 2) * self.towhead_mass * (self.radius ** 2))
        self.Iy = 200
        self.Iz = self.Iy

        # Complete mass matrix for rigid body
        self.mass_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom))
        for i in range(3):
            self.mass_matrix[i][i] = self.mass
        self.mass_matrix[3][3] = self.Ix
        self.mass_matrix[4][4] = self.Iy
        self.mass_matrix[5][5] = self.Iz

        for i in range(self.resonant_mode_count):
            for j in range(self.resonant_mode_count):
                self.mass_matrix[6+i][6+j] = 2


        self.damping_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom))
        self.stiffness_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom))

    def substitutions(self):
        return {'z_s': self.depth, 'mode_count': self.degrees_of_freedom, 'seafloor_depth': self.seafloor_depth,
                'mass_matrix': self.mass_matrix, 'damping_matrix': self.damping_matrix,
                'stiffness_matrix': self.stiffness_matrix
                }

    def geometry(self):
        geometry = pygmsh.opencascade.Geometry()
        geometry.add_cylinder(x0=[-self.length / 2, 0, self.depth], axis=[self.length, 0, 0],
                              radius=self.radius, angle=2 * np.pi, char_length=1)
        return geometry


# Design variables
print('Case: ', str(1))
design_variables = np.asarray([-1, 0.274, 0.01, 10, 0, 0, 110])
tube = FlexibleTube(design_variables)
tube_substitutions = tube.substitutions()
tube_geometry = tube.geometry()
create_case_directory(output_directory, 1)
create_case_files(common_file_directory, tube_substitutions)
create_mesh_file(tube_geometry, device_name, gmsh_exe_location, mesh_refinement_factor)
submerged_mesh(device_name)
#run_wamit(run_wamit_command, flexible_bool=True)
print('Done.')
