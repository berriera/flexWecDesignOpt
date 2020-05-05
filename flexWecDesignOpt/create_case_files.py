def create_case_files(common_bem_file_folder, copy_folder, design_vars):
    """Copies original input files from their common directory, and changes them in the process using variable
    substitution.

    Args:
        common_bem_file_folder (str):
        copy_folder (str):
        design_vars (one-dimensional array):

    Returns:
        None
    """

    def change_case_file(text_file, design_var):
        with open(text_file, 'r') as g:
            line_list = []
            line_number = 0
            lines_text = g.readlines()
            for line in lines_text:
                if line.find('?') != -1:
                    index = 0
                    replace_keys = []
                    for character in line:
                        if character == '?':
                            replace_keys.append(index)
                        index += 1
                    replace_count = int(len(replace_keys) / 2)
                    for replacement in range(0, replace_count):
                        string_beginning = replace_keys[replacement]
                        string_end = replace_keys[replacement + 1]
                        replacement_variable_number = int(line[string_beginning + 1:string_end])
                        print(replacement_variable_number)
                        old_string = line[string_beginning:string_end + 1]
                        print('old_string')
                        print(old_string)
                        replacement_string = str(design_var[replacement_variable_number - 1])
                        print('replacement_string')
                        print(replacement_string)
                        line = line.replace(old_string, replacement_string)
                line_list.append(line)
                line_number += 1
            return line_list

    import shutil
    import os
    files = os.listdir(common_bem_file_folder)
    for bem_input_file in files:
        file = common_bem_file_folder + '/' + bem_input_file
        new_file_text = change_case_file(file, design_vars)
        shutil.copy(file, copy_folder)
        change_file = copy_folder + '/' + bem_input_file
        with open(change_file, 'w') as f:
            for new_string in new_file_text:
                f.write(new_string)
    return
