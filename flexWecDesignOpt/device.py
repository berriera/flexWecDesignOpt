class Device(object):

    # Simulation information
    name = 'device'
    boundary_element_method = 'wamit'
    verbose = True

    # File information
    input_file_location = ''
    copy_file_list = ['']
    edit_file_list = ['']
    output_file_location = ''

    # Boundary element method file information
    gmsh_location = ''
    wamit_location = r'C:\WAMITv7\wamit'
    defmod_location = None

    # Objective weight information
    weights = []

    def __init__(self, design_variables=None):
        self.substitution_dict = {}
        self.mesh_geometry = None
        self.design_variables = design_variables

        self.folder_name = self.name
        self.case_directory = self.output_file_location

        self.objectives = []

    def parse_input(self):
        """Opens and reads in specified parameters in the user created .yaml file

        Args:
            input_file_location (str): location of input.yaml file

        Returns
            file_names (dict): dictionary of keys: items from input.yaml file
        """
        import yaml

        with open(self.input_file_location, 'r') as f:
            file_names = yaml.safe_load(f)
        return file_names

    def create_case_directory(self):
        import os

        # Count how many designs have already been evaluated in the output folder
        output_folder_file_list = os.listdir(self.output_file_location)
        folder_count = 1
        for file in output_folder_file_list:
            if file.startswith(self.folder_name):
                folder_count += 1
        self.case_directory = os.path.join(self.output_file_location, self.folder_name + '_' + str(folder_count))

        # Creates the case output folder if it does not already exist
        assert not os.path.exists(self.case_directory), "Design output folder already exists"  # TODO: remove?
        os.makedirs(self.case_directory)
        os.chdir(self.case_directory)

    def create_case_files(self, substitution_dict=None):
        import shutil
        import os

        import write_input_files

        if self.verbose:
            print('\tWriting input files...')

        # Make list of all files in folder to be copied or edited
        file_copy_list = Device.copy_file_list + Device.edit_file_list

        # Copy each input file into the
        for input_file in file_copy_list:
            full_file_name = os.path.join(self.input_file_location, input_file)
            if input_file in file_copy_list:
                shutil.copy(full_file_name, self.case_directory)

                # Edit file with substitution variables
                if input_file in Device.edit_file_list:
                    new_file_text = write_input_files.change_case_file(full_file_name, substitution_dict)
                    shutil.copy(full_file_name, self.case_directory)
                    change_file = os.path.join(self.case_directory, input_file)
                    with open(change_file, 'w') as f:
                        for new_line in new_file_text:
                            f.write(new_line)

        if self.verbose:
            print('\tWriting complete.')

        return

    def mesh(self):
        import os
        import pygmsh

        if self.verbose:
            print('\tMeshing...')

        assert os.path.isfile(self.gmsh_location) or os.path.isfile(self.gmsh_location + '.exe'), \
            'GMSH.exe file in location {} not found'.format(self.gmsh_location)

        mesh = pygmsh.generate_mesh(self.mesh_geometry,
                                    dim=2,
                                    remove_lower_dim_cells=True,
                                    gmsh_path=self.gmsh_location,
                                    mesh_file_type='stl',
                                    verbose=False)

        mesh_file_name = self.name + '.stl'
        mesh.write(mesh_file_name)
        mesh_stl_location = os.path.join(self.case_directory, mesh_file_name)

        if self.boundary_element_method.lower() == 'wamit':
            self._write_gdf(mesh_stl_location)

        if self.verbose:
            print('\tMeshing complete.')
        return

    def _write_gdf(self, mesh_stl_location):
        import numpy as np
        import meshmagick.mesh_clipper
        import meshmagick.mmio
        import meshmagick.mesh

        vertices, faces = meshmagick.mmio.load_STL(mesh_stl_location)
        if np.any(vertices[:, 2] >= 0.0):
            mesh_all = meshmagick.mesh.Mesh(vertices, faces)
            mesh_clip_object = meshmagick.mesh_clipper.MeshClipper(source_mesh=mesh_all)
            mesh_clipped_below_waterline = mesh_clip_object.clipped_mesh
            vertices = mesh_clipped_below_waterline.vertices
            faces = mesh_clipped_below_waterline.faces
        meshmagick.mmio.write_GDF(self.name + '.gdf', vertices, faces)
        # TODO: healnormals function on mesh_geometry
        # TODO: user specified mesh_geometry symmetry arguments for quicker write and meshing times for WAMIT using a
        #  symmetry list device argument

    def run_wamit(self):
        import os
        import subprocess

        os.chdir(self.case_directory)

        if self.verbose:
            print('\tRunning analysis...')

        run_executable_list = [Device.wamit_location]
        if Device.defmod_location is not None:
            run_executable_list.append(Device.defmod_location)
            run_executable_list.append(Device.wamit_location)

        for analysis_executable in run_executable_list:
            if self.verbose:
                print('\tRunning ' + os.path.basename(analysis_executable))
            assert os.path.isfile(analysis_executable), 'Analysis .exe file in location not found at {}'.format(
                analysis_executable)
            subprocess.run([analysis_executable])

        if self.verbose:
            print('\tAnalysis completed.')
        return

    def run_nemoh(self):
        pass

    def read_wamit(self):
        import matlab.engine
        import scipy.io

        # Read WAMIT output file using matlab API
        output_file_name = self.case_directory + Device.name + '.out'
        matlab_engine = matlab.engine.start_matlab()
        matlab_engine.read_WAMIT_output(output_file_name)

        # Load .mat file of generated results and send them to a dictionary
        results_dict = {}
        scipy.io.loadmat(self.case_directory + 'WAMIT_results.mat', mdict=results_dict, squeeze_me=True)
        results_dict = results_dict['hydro']

        return results_dict

    def delete_case_files(self):  # TODO: create function to delete input_files and output files of given names
        return

    def evaluate_objective(self):
        pass
