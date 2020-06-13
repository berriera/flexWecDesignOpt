def create_mesh_file(geometry, device_name, gmsh_exe_location, mesh_refinement_factor=0.5):
    """This function creates the initial mesh file in .stl format using the software GMSH

    Args:
        geometry (geometry object): generated from the device.geometry class from a set of design variables
        device_name (str): name of the device
        gmsh_exe_location (str): file location of the installed GMSH software (should end in 'gmsh')
        mesh_refinement_factor (float): specifies how refined the generated mesh is. Lower means more refined.
                                            Default is 0.5

    Returns:
        None

    """
    import pygmsh

    print('\tMeshing...')

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
    """This function clips the .stl mesh below the waterline then writes the new .gdf mesh file to the current
    output folder.

    Args:
        device_name (str): name of the device

    Returns:
        None

    """
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
