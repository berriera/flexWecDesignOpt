import pygmsh
import os
import numpy as np
from file_mgmt import create_case_directory
from mesh import create_mesh_file
from mesh import submerged_mesh

gmsh_exe_location = 'C:/Users/13365/Documents/gmsh-4.5.6-Windows64/gmsh-4.5.6-Windows64/gmsh'
output_directory = os.path.abspath(os.path.join('examples', 'output'))
device_name = 'S3'
mesh_refinement_factor = 0.45

# Design variables
depth = -1
radius = 0.274
thickness = 0.01
length = 6
end_cap_thickness = thickness

geometry = pygmsh.opencascade.Geometry()
outer_tube = geometry.add_cylinder(x0=[0, 0, depth], axis=[length, 0, 0], radius=radius,
                                   angle=2 * np.pi, char_length=1)
inner_tube = geometry.add_cylinder(x0=[end_cap_thickness, 0, depth], radius=radius - thickness,
                                   axis=[length - 2 * end_cap_thickness, 0, 0], angle=2 * np.pi, char_length=1)
geometry.boolean_difference([outer_tube], [inner_tube])

create_case_directory(output_directory, 3)
create_mesh_file(geometry, device_name, gmsh_exe_location, mesh_refinement_factor)
submerged_mesh(device_name)
