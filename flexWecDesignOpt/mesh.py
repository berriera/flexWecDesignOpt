def create_mesh_file(geometry, device_name, gmsh_exe_location, mesh_refinement_factor=0.5):
    import pygmsh  # TODO: documentation

    meshing_arguments = ['-clscale', str(mesh_refinement_factor),  # set mesh element size factor
                         '-clcurv', str(360 / 50),  # computes mesh element size from curvature
                         '-setnumber', 'Mesh.SubdivisionAlgorithm', '1',
                         # subdivision algorithm 1 means all quadrangles
                         '-setnumber', 'Mesh.RecombineAll', '1']  # applies recombination algorithm to all surfaces

    mesh = pygmsh.generate_mesh(geometry,
                                dim=2,
                                extra_gmsh_arguments=meshing_arguments,
                                remove_lower_dim_cells=True,
                                gmsh_path=gmsh_exe_location,
                                geo_filename=device_name + '.geo',
                                msh_filename=device_name + '.stl',
                                mesh_file_type='stl',
                                verbose=False)


def submerged_mesh(device_name):
    # Clip the .stl mesh below the waterline then write the new .gdf mesh file
    import meshmagick.mesh_clipper
    import meshmagick.mmio
    import meshmagick.mesh
    vertices, faces = meshmagick.mmio.load_STL(device_name + '.stl')
    mesh_all = meshmagick.mesh.Mesh(vertices, faces)  # TODO: healnormals function on mesh
    mesh_clip_object = meshmagick.mesh_clipper.MeshClipper(source_mesh=mesh_all)
    mesh_clipped_below_waterline = mesh_clip_object.clipped_mesh
    meshmagick.mmio.write_GDF(device_name + '.gdf', mesh_clipped_below_waterline.vertices,
                              mesh_clipped_below_waterline.faces)
    # TODO: user specified mesh symmetry arguments for quicker write and meshing times,
    #  need to change meshmagick for this
