def create_mesh_file(geometry, device_name, gmsh_exe_location, verbosity=True):
    # TODO: set default device_name = 'mesh_file'
    """This function creates the initial mesh file in .stl format using the software GMSH

    Args:
        geometry (geometry object): generated from the device.geometry class from a set of design variables
        device_name (str): name of the device
        gmsh_exe_location (str): file location of the installed GMSH software (should end in 'gmsh')
        verbosity (bool): controls printing verbosity for writing mesh file

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

    mesh_file_name = device_name + '.stl'
    mesh.write(mesh_file_name)
    _trim_mesh_below_waterline(device_name=device_name)

    if verbosity:
        print('\tDone.')


def _trim_mesh_below_waterline(device_name):
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

    # TODO: make writing mesh file BEM software dependent (WAMIT, NEMOH, Capytaine)
    # TODO: user specified mesh symmetry arguments for quicker write and meshing times,
    #  need to change meshmagick for this
    return
