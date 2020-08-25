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

#common_file_directory = os.path.abspath(common_file_directory)
#output_directory = os.path.abspath(output_directory)

# Design variables
depth = -1
radius = 0.274
thickness = 0.01
length = 10
end_cap_thickness = thickness

geometry = pygmsh.opencascade.Geometry()
outer_tube = geometry.add_cylinder(x0=[-length / 2, 0, depth], axis=[length, 0, 0], radius=radius,
                                   angle=2 * np.pi, char_length=1)

create_case_directory(output_directory, 1)
create_mesh_file(geometry, device_name, gmsh_exe_location, mesh_refinement_factor)
submerged_mesh(device_name)
run_wamit(run_wamit_command, flexible_bool=True)
