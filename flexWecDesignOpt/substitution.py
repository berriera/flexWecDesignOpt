def substitute_variables_in_line(line_text, variables, file_name='', line_number=1):  # TODO: documentation
    """Alters a single line in a text file with included ?var? statements

    Args:
        line_text (str): A single line in a text file to be altered
        variables (one-dimensional array or dict): Design variable array or dict of variable substitutions
        file_name (str): Name of the file. Used for error handling to point user to substitution mistake
        line_number (int): Line number of file_name. Used for error handling to point user to variable substitution
                            mistake

    Returns:
        line_text (str) : The substituted line text for the file

    """
    import numpy as np
    error_list = []  # TODO: only repeat key error once (fix error_list var to keep appending to existing error_list
    while line_text.find('?') != -1:
        substitution_indices = [position for position, character in enumerate(line_text) if character == '?']
        if len(substitution_indices) % 2 != 0.0:
            raise ValueError("Check line number " + str(line_number) + " in file  " + file_name +
                             " for proper substitution formatting")
        string_substitution = line_text[substitution_indices[0]:substitution_indices[1] + 1]
        if variables is None:
            return
        elif isinstance(variables, np.ndarray):
            # TODO: enforce try else statements here for checking if var is in dict or if array is out of bounds
            replacement_variable_number = int(line_text[substitution_indices[0] + 1: substitution_indices[1]])
            replacement_string = str(variables[replacement_variable_number - 1])
        elif isinstance(variables, dict):
            try:
                replacement_variable_key = line_text[substitution_indices[0] + 1: substitution_indices[1]]
                replacement_string = str(variables[replacement_variable_key])
            except KeyError:
                error_message = "\tError: Check substitution method for missing key " + replacement_variable_key + \
                                ".\n\tIt is found on line number " + str(line_number) + " in input file " + file_name
                if error_message not in error_list:
                    error_list.append(error_message)
                    print(error_message)
        else:
            raise TypeError("Returned object type of device substitution should be an array or a dictionary.")
        try:
            line_text = line_text.replace(string_substitution, replacement_string)
        except UnboundLocalError:
            print("Could not replace variable. Check input file for corrections.")
            break
    return line_text


def change_case_file(text_file, design_var):
    """Alters a single input file in a directory

    Args:
        text_file (str): Name of the text file to be read and potentially altered
        design_var (one-dimensional array or dict): Design variable array or dict of variable substitutions

    Returns:
        line_list (str list): list of altered lines
    """
    with open(text_file, 'r') as g:
        line_list = []
        line_number = 1
        lines_text = g.readlines()
        for line in lines_text:
            if line.find('?') != -1:
                line = substitute_variables_in_line(line, design_var, text_file, line_number)
            line_list.append(line)
            line_number += 1
        return line_list


def create_case_files(common_file_directory, substitution_array):
    """Copies original input files from their common directory, and changes them in the process using variable
    substitution.

    Args:
        common_file_directory (str): Input file folder location
        case_output_folder (str): Output file folder location
        substitution_array (one-dimensional array or dict): Design variable array or dict of variable substitutions

    Returns:
        None
    """

    import shutil
    import os
    extension_do_not_copy_list = ['.yaml', 'csv', '.py']
    files = os.listdir(common_file_directory)
    for bem_input_file in files:
        file = os.path.join(common_file_directory, bem_input_file)
        extension = os.path.splitext(file)[1]
        if extension not in extension_do_not_copy_list:
            new_file_text = change_case_file(file, substitution_array)
            shutil.copy(file, os.getcwd())
            change_file = os.path.join(os.getcwd(), bem_input_file)
            with open(change_file, 'w') as f:
                for new_string in new_file_text:
                    f.write(new_string)
    return
