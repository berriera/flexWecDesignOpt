def create_case_directory(case_number, output_directory):
    """Creates a directory to hold each design's input and generated output files.

    Args:
        case_number (int):
        output_directory (str):

    Returns:
        str
    """
    import os
    path = os.path.join(output_directory, 'case_' + str(case_number))
    os.makedirs(path)
    return path


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
