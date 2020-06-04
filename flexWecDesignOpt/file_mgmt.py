def create_case_directory(output_directory='', case_number=1):
    """Creates a directory to hold each design's input and generated output files.

    Args:
        case_number (int):
        output_directory (str):

    Returns:
        str
    """
    import os
    if case_number >= 0:
        file_path_list = os.listdir(output_directory)
        folder_count = len(file_path_list)
        directory_path = os.path.join(output_directory, 'case_' + str(case_number))
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)
        else:
            case_number = folder_count + 1
            print('\tCase output folder already exists. Renaming to case_' + str(case_number))
            directory_path = os.path.join(output_directory, 'case_' + str(case_number))
            os.makedirs(directory_path)
    return directory_path


def current_case_directory():
    return create_case_directory(case_number=-1)


def parse_input(input_file_location):
    """Opens and reads in specified parameters in the user created .yaml file

    Args:
        input_file_location (str):

    Returns
        str
    """
    import yaml
    with open(input_file_location, 'r') as f:
        file_names = yaml.safe_load(f)
    return file_names
