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

    """  # TODO: comment code
    import numpy as np

    # Search for ? symbols in a given line of an input file
    while line_text.find('?') != -1:
        substitution_indices = [position for position, character in enumerate(line_text) if character == '?']

        # Check that each substitution variable has two surrounding question marks
        if len(substitution_indices) % 2 != 0.0:
            raise ValueError("Check line number " + str(line_number) + " in file  " + file_name +
                             " for proper substitution formatting")
        string_substitution = line_text[substitution_indices[0]:substitution_indices[1] + 1]
        if variables is None:  # TODO: remove line after changing substitution file list
            return
        elif isinstance(variables, np.ndarray):  # TODO: remove array substitution functionality
            # TODO: check if var is in dict or if array is out of bounds
            replacement_variable_number = int(line_text[substitution_indices[0] + 1: substitution_indices[1]])
            replacement_string = str(variables[replacement_variable_number - 1])
        elif isinstance(variables, dict):
            try:  # TODO: replace try catch with assertion errors
                replacement_variable_key = line_text[substitution_indices[0] + 1: substitution_indices[1]]
                replacement_variable_value = variables[replacement_variable_key]

                # Substitute real number in line
                if isinstance(replacement_variable_value, int) or isinstance(replacement_variable_value, float):
                    replacement_string = str(variables[replacement_variable_key])

                # Substitute 1d or 2d np array in input file
                elif isinstance(replacement_variable_value, np.ndarray):
                    if replacement_variable_value.ndim == 1:
                        array_string = np.array2string(replacement_variable_value, max_line_width=10000)
                        replacement_string = array_string[1:-1]

                    elif replacement_variable_value.ndim == 2:
                        replacement_string = ""
                        row_count = replacement_variable_value.shape[0]
                        for row in range(row_count):
                            matrix_row = replacement_variable_value[row]
                            matrix_row_string = np.array2string(matrix_row, max_line_width=10000)
                            replacement_string = replacement_string + matrix_row_string[1:-1]
                            if row != row_count - 1:
                                replacement_string = replacement_string + '\n'
            except KeyError:
                error_message = "\tError: Check substitution method for missing key " + replacement_variable_key + \
                                ".\n\tIt is found on line number " + str(line_number) + " in input file " + file_name
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
