import pygmsh
import os
import numpy as np
from file_mgmt import create_case_directory
from mesh import create_mesh_file
from mesh import submerged_mesh
from analysis import run_wamit

run_wamit_command = 'C:\WAMITv7\wamit'
gmsh_exe_location = 'C:/Users/13365/Documents/gmsh-4.5.6-Windows64/gmsh-4.5.6-Windows64/gmsh'
common_file_directory = os.path.abspath(os.path.join('examples', 'flexible_tube'))
output_directory = os.path.abspath(os.path.join('examples', 'output'))
device_name = 'tube'
mesh_refinement_factor = 0.40


class FlexibleTube(object):

    def __init__(self, design_vars):
        self.device_name = 'tube'
        self.depth = design_vars[0]
        self.radius = design_vars[1]
        self.thickness = design_vars[2]
        self.length = design_vars[3]
        self.distensibility = design_vars[4]
        self.eta = design_vars[5]

    def substitutions(self):
        return {'z_s': self.depth, 'r_s': self.radius, 't': self.thickness, 'L': self.length
                }

    def geometry(self):
        geometry = pygmsh.opencascade.Geometry()
        geometry.add_cylinder(x0=[-self.length / 2, 0, self.depth], axis=[self.length, 0, 0],
                              radius=self.radius, angle=2 * np.pi, char_length=1)
        return geometry


# Design variables
print('Case: ', str(1))
design_variables = np.asarray([-1, 0.274, 0.01, 10, 0, 0])
tube = FlexibleTube(design_variables)
tube_substitutions = tube.substitutions()
tube_geometry = tube.geometry()
create_case_directory(output_directory, 1)
create_mesh_file(tube_geometry, device_name, gmsh_exe_location, mesh_refinement_factor)
submerged_mesh(device_name)
run_wamit(run_wamit_command, flexible_bool=True)
print('Done.')