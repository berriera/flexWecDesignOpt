import argparse
import sys
import numpy as np

# sys.path.append(os.path.dirname(os.path.realpath(__file__)))
# import flexWecDesignOpt

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
                    required=True,
                    help='file path to input.yaml file'
                    )
parser.add_argument('-r', '--run',
                    action='store_true',
                    default=False,
                    help='run boundary element method software on created input files'
                    )
args = parser.parse_args(sys.argv[1:])


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


def create_case_directory(case_number, output_directory):
    """Creates a directory to hold each design's input and generated output files.

    Args:
        case_number (int):
        output_directory (str):

    Returns:
        str
    """
    import os
    path = output_directory + '/case_' + str(case_number)
    os.makedirs(path)
    return path


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


def run_wamit(output_folder, bem_command):
    """
    This function runs the boundary element solver command in the current case directory.

    Args:
        output_folder(str):
        bem_command (str):

    Returns:
        None
    """
    import os
    import subprocess
    os.chdir(output_folder)
    subprocess.run([bem_command])
    return


def read_output():
    """

    :return:
    """
    return


def main():
    input_file_names = parse_input(args.input)
    cases_file = input_file_names['cases_file']
    design_data = np.genfromtxt(cases_file, delimiter=',')
    design_count = design_data.shape[0]
    case_output_folder = input_file_names['output_directory']

    for case in range(design_count):
        print(case)
        design_variables = design_data[case, :]
        print(design_variables)
        case_output_folder = create_case_directory(case + 1, input_file_names['output_directory'])
        create_case_files(input_file_names['common_file_directory'], case_output_folder, design_variables)
        if args.run:
            run_wamit(case_output_folder, input_file_names['run_wamit_directory'])
            read_output()


if __name__ == "__main__":
    main()
