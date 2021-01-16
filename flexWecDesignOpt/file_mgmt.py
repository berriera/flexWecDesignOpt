def create_case_directory(output_directory='', folder_name='design_case_', analysis_type=''):
    """Creates a directory to hold each design's input files and generated output files.

    Args:
        output_directory (str): folder location for all generated results
        folder_name (str): folder naming convention
        analysis_type (str): keyword for analysis software

    Returns:
        directory_path (str): location of the created case directory
    """
    import os

    # Count how many designs have already been evaluated in the output folder
    file_path_list = os.listdir(output_directory)
    folder_count = 1
    for file in file_path_list:
        if file.startswith(folder_name):
            folder_count += 1
    case_folder_path = os.path.join(output_directory, folder_name + str(folder_count))

    # Creates the case output folder if it does not already exist
    assert not os.path.exists(case_folder_path), "Design output folder already exists"  # TODO: remove?
    os.makedirs(case_folder_path)
    subdirectories = file_information(analysis_type)['subdirectories']
    for directory in subdirectories:
        subdirectory_path = os.path.join(case_folder_path, directory)
        os.makedirs(subdirectory_path)
    os.chdir(case_folder_path)
    return case_folder_path


def delete_case_files():  # TODO: create function to delete input_files and output files that uses file extensions
    return


def parse_input(input_file_location):
    """Opens and reads in specified parameters in the user created .yaml file

    Args:
        input_file_location (str): location of input.yaml file

    Returns
        file_names (dict): dictionary of keys: items from input.yaml file
    """
    import yaml
    with open(input_file_location, 'r') as f:
        file_names = yaml.safe_load(f)
    return file_names


def file_information(analysis_software_key):
    """Returns file information for each analysis type.
    'subdirectory' key contains a list of required subdirectory names in a design's folder
    'input_extensions' key contains a list of input file extensions

    """
    analysis_software_key = analysis_software_key.lower()

    file_info = {
        'wamit': {'mesh_subdirectory': None,
                  'bem_subdirectory': None,
                  'subdirectories': [],
                  'input_extensions': [],
                  },
        'nemoh': {'mesh_subdirectory': None,
                  'bem_subdirectory': None,
                  'subdirectories': [],
                  'input_extensions': []
                  },
        'wec-sim': {'mesh_subdirectory': 'geometry',
                    'bem_subdirectory': 'hydroData',
                    'subdirectories': ['hydroData', 'geometry'],
                    'input_extensions': []
                    }
    }

    # TODO: update for NEMOH
    # TODO: figure out how this should work for using both WAMIT and WEC-Sim
    return file_info[analysis_software_key]
