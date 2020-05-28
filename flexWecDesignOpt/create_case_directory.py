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
