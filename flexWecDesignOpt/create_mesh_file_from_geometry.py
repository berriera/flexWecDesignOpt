def create_mesh_file_from_geometry(geometry, device_name, case_output_folder,
                                   gmsh_exe_location, mesh_refinement_factor=0.5):
    import os
    import pygmsh
    import meshmagick.mesh_clipper
    import meshmagick.mmio
    import meshmagick.mesh

    os.chdir(case_output_folder)

    meshing_arguments = ['-clscale', str(mesh_refinement_factor),  # set mesh element size factor
               '-clcurv', str(360 / 50),  # computes mesh element size from curvature
               '-setnumber', 'Mesh.SubdivisionAlgorithm', '1',  # subdivision algorithm 1 means all quadrangles
                '-setnumber', 'Mesh.RecombineAll', '1'] # applies recombination algorithm to all surfaces

    mesh = pygmsh.generate_mesh(geometry,
                                dim=2,
                                extra_gmsh_arguments=meshing_arguments,
                                remove_lower_dim_cells=True,
                                gmsh_path=gmsh_exe_location,
                                geo_filename=device_name + '.geo',
                                msh_filename=device_name + '.stl',
                                mesh_file_type='stl',
                                verbose=False)

    # Clip the .stl mesh below the waterline then write the new .gdf mesh file
    # TODO: check normals
    # TODO: make all faces end at zero
    vertices, faces = meshmagick.mmio.load_STL(device_name + '.stl')
    mesh_all = meshmagick.mesh.Mesh(vertices, faces)
    mesh_clip = meshmagick.mesh_clipper.MeshClipper(source_mesh=mesh_all)
    mesh_submerged = mesh_clip.lower_mesh
    meshmagick.mmio.write_GDF(device_name + '.gdf', mesh_submerged.vertices, mesh_submerged.faces)
