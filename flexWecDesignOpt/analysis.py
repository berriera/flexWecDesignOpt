def run_analysis(executable_locations, analysis_type, verbosity=True):
    """ # TODO: update documentation
    This function runs the boundary element solver command in the current case directory.

    Args:
        executable_locations (list of raw str): location of boundary element method executable files
        verbosity (bool): specify if print statements are outputted to screen
    Returns:
        None
    """
    import os
    import subprocess

    import file_mgmt

    if verbosity:
        print('\tRunning analysis...')

    current_folder = os.getcwd()
    relative_bem_output_path = file_mgmt.file_information(analysis_type)['bem_subdirectory']
    if relative_bem_output_path is not None:
        absolute_bem_folder_path = os.path.join(current_folder, relative_bem_output_path)
        os.chdir(absolute_bem_folder_path)

    for analysis_executable in executable_locations:
        if verbosity:
            print('\tRunning ' + os.path.basename(analysis_executable))
        assert os.path.isfile(analysis_executable), 'Analysis .exe file in location not found'
        subprocess.run([analysis_executable])

    if relative_bem_output_path is not None:
        os.chdir(current_folder)

    if verbosity:
        print('\tDone.')
    return


def boundary_condition_frequency_solver(function_name, eigenvalue_count=1, h=1.0, freq_limit=50, reltol=1e-5, eps=1e-3,
                                        extra_args=()):
    """
    This function finds a user-specified number of eigenvalue solutions to a device's characteristic equation.

    Args:
        function_name: the name of the system to be iteratively solved
        eigenvalue_count (int): number of unique roots to f(x) = 0 be found
        h (float): increment used to guess new roots; default is 1.0
        freq_limit (float): maximum allowable frequency to be found; used as a break in case of unsuccessful root
                            finding or to only search a certain space; default is 100
        reltol (float): relative tolerance used to find if a new found root is unique; default is 1e-5
        eps(float): used to as a minimum limit on found roots to avoid solutions too close to 0; default is 1e-3
        extra_args (tuple): used to pass extra constant values to fsolve

    Returns:
        nonzero_eigenvalues (numpy array): list of non-zero roots to the equation; shape is 1 row by eigenvalue_count
                                            columns or smaller if less roots in the specified range were found
    """

    import numpy as np
    import scipy.optimize

    x0 = reltol  # used to avoid division by zero errors in solving certain problems with fsolve
    found_eigenvalue_count = 0
    eigenvalues = np.zeros(eigenvalue_count, )
    while found_eigenvalue_count < eigenvalue_count and x0 < freq_limit:
        root = scipy.optimize.fsolve(function_name, x0, extra_args)[0]
        if root > eps and not np.any(abs(eigenvalues - root) / root < reltol):
            eigenvalues[found_eigenvalue_count] = root
            found_eigenvalue_count = found_eigenvalue_count + 1
        x0 = x0 + h
    nonzero_eigenvalues = eigenvalues[eigenvalues != 0]
    return nonzero_eigenvalues


# TODO: rewrite to travel along function_name and make brackets so there's no need for h

def evaluate_device(device_class, design_variables):
    # External import
    import os

    # Internal imports
    import write_input_files
    import mesh
    import analysis  # TODO: fix this to not import own file (move to new file entirely called evaluate.py)
    import file_mgmt
    import output

    # Create a device object
    device_object = device_class(design_variables)
    input_file_substitutions = device_object.substitutions
    mesh_geometry = device_object.mesh

    # Create output folder for device evaluation
    working_folder = file_mgmt.create_case_directory(output_directory=device_class.output_file_directory,
                                                     analysis_type=device_class.analysis_type)

    # Create input_files files in directory and then run specified analysis executables in that directory
    write_input_files.create_case_files(device_object.input_file_directory, input_file_substitutions)
    mesh.create_mesh_file(device_object=device_object,
                          verbosity=device_object.verbosity)
    analysis.run_analysis(executable_locations=device_class.exe_locations,
                          analysis_type=device_class.analysis_type,
                          verbosity=device_object.verbosity)

    # Output logger
    output.write_dict_to_text_file(input_file_substitutions)
    os.chdir(device_object.input_file_directory)
    objective_function = device_object.objective(working_folder)
    return objective_function
