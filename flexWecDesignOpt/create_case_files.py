def string_change(line_change, variables, file_name='', line_number=1):  # TODO: documentation
    while line_change.find('?') != -1:
        substitution_indices = [position for position, character in enumerate(line_change) if character == '?']
        if len(substitution_indices) % 2 != 0.0:
            raise ValueError("Check line number" + str(line_number) + "in file" + file_name +
                             "for proper substitution formatting")
        string_substitution = line_change[substitution_indices[0]:substitution_indices[1] + 1]
        if isinstance(variables, list):
            replacement_variable_number = int(line_change[substitution_indices[0] + 1: substitution_indices[1]])
            replacement_string = str(variables[replacement_variable_number - 1])
        elif isinstance(variables, dict):
            replacement_variable_key = line_change[substitution_indices[0] + 1: substitution_indices[1]]
            replacement_string = line_change[replacement_variable_key]
        else:
            raise TypeError("Returned object type of device substitution should be an array or a dictionary.")
        line_change = line_change.replace(string_substitution, replacement_string)
    return line_change


def change_case_file(text_file, design_var):  # TODO: documentation
    with open(text_file, 'r') as g:
        line_list = []
        line_number = 0
        lines_text = g.readlines()
        for line in lines_text:
            if line.find('?') != -1:
                line = string_change(line, design_var)
            line_list.append(line)
            line_number += 1
        return line_list


def create_case_files(common_file_directory, case_output_folder, substitution_array):
    """Copies original input files from their common directory, and changes them in the process using variable
    substitution.

    Args:
        common_file_directory (str):
        case_output_folder (str):
        substitution_array (one-dimensional array):

    Returns:
        None
    """

    import shutil
    import os
    files = os.listdir(common_file_directory)
    for bem_input_file in files:
        file = common_file_directory + '/' + bem_input_file  # TODO: exclude .csv files
        new_file_text = change_case_file(file, substitution_array)
        shutil.copy(file, case_output_folder)
        change_file = case_output_folder + '/' + bem_input_file  # TODO: join these strings in a non-Windows specific way
        with open(change_file, 'w') as f:
            for new_string in new_file_text:
                f.write(new_string)
    return
