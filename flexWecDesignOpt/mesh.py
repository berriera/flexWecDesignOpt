def create_mesh_file(geometry, device_name, gmsh_exe_location, verbosity=True):
    # TODO: set default device_name = 'mesh_file'
    # TODO: make submerged mesh internal
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

    if verbosity:
        print('\tMeshing...')

    mesh = pygmsh.generate_mesh(geometry,
                                dim=2,
                                remove_lower_dim_cells=True,
                                gmsh_path=gmsh_exe_location,
                                mesh_file_type='stl',
                                verbose=False)
    mesh.write(device_name + '.stl')


def submerged_mesh(device_name):
    # TODO: change default device name, adjust analysis.py accordingly
    """This function clips the .stl mesh below the waterline then writes the new .gdf mesh file to the current
    output folder.

    Args:
        device_name (str): name of the device

    Returns:
        None

    """
    import numpy as np
    import meshmagick.mesh_clipper
    import meshmagick.mmio
    import meshmagick.mesh
    vertices, faces = meshmagick.mmio.load_STL(device_name + '.stl')
    if np.any(vertices[:, 2] >= 0):
        mesh_all = meshmagick.mesh.Mesh(vertices, faces)  # TODO: healnormals function on mesh
        mesh_clip_object = meshmagick.mesh_clipper.MeshClipper(source_mesh=mesh_all)
        mesh_clipped_below_waterline = mesh_clip_object.clipped_mesh
        vertices = mesh_clipped_below_waterline.vertices
        faces = mesh_clipped_below_waterline.faces
    meshmagick.mmio.write_GDF(device_name + '.gdf', vertices, faces)
    return vertices
    # TODO: user specified mesh symmetry arguments for quicker write and meshing times,
    #  need to change meshmagick for this
