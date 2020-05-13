import pygmsh
import numpy as np
import meshmagick.mesh_clipper
import meshmagick.mmio
import meshmagick.mesh

geom = pygmsh.opencascade.Geometry()

fbname = 'cylinder_test'
ofst = 0.1
geom.add_cylinder(x0=[0.0, 0.0, 0+ofst],
                         axis=[0.0, 0.0, -2],
                         radius=1,
                         angle=2 * np.pi,
                         char_length=1)

mshRefFactor = 0.25

mshArgs = ['-clscale', str(mshRefFactor),                   # set mesh element size factor
           '-clcurv', str(360/50),                          # computes mesh element size from curvature
           '-setnumber', 'Mesh.SubdivisionAlgorithm', '1', # subdivision algorithm 1 means all quadrangles
           '-setnumber', 'Mesh.RecombineAll', '1']          # applies recombination algorithm to all surfaces

mesh = pygmsh.generate_mesh(geom,
                            dim=2,
                            extra_gmsh_arguments=mshArgs,
                            remove_lower_dim_cells=True,
                            gmsh_path='C:/Users/13365/Documents/gmsh-4.5.6-Windows64/gmsh-4.5.6-Windows64/gmsh',
                            geo_filename=fbname + '.geo',
                            msh_filename=fbname + '.stl',
                            mesh_file_type='stl')

vertices, faces = meshmagick.mmio.load_STL('cylinder_test.stl')
mesh_all = meshmagick.mesh.Mesh(vertices, faces)
mesh_clip = meshmagick.mesh_clipper.MeshClipper(source_mesh=mesh_all)
mesh_submerged = mesh_clip.lower_mesh
meshmagick.mmio.write_GDF('cylinder_test_full.gdf', mesh_all.vertices, mesh_all.faces)
meshmagick.mmio.write_GDF('cylinder_wamit_submerged.gdf', mesh_submerged.vertices, mesh_submerged.faces)
