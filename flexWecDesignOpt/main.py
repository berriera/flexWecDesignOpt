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
                    default=True,
                    type=bool,
                    help='run boundary element method software on created input files'
                    )
args = parser.parse_args(sys.argv[1:])


def parse_input(input_file_location):
    """

    :return:
    """
    import yaml
    with open(input_file_location, 'r') as f:
        file_names = yaml.safe_load(f)
    return file_names


def create_case_directory(case_number, output_directory):
    """

    :param case_number:
    :param output_directory:
    :return:
    """
    import os
    path = output_directory + '/case_' + str(case_number)
    os.makedirs(path)
    return path


def create_case_files(common_bem_file_folder, copy_folder, design_vars):
    """

    :param common_bem_file_folder:
    :param copy_folder:
    :return:
    """

    def change_case_file(text_file, design_var):
        with open(text_file, 'r') as g:
            line_list = []
            for line in text_file:
                line_text = g.readline(line)
                if line_text.find('?') != -1:
                    index = 0
                    replace_keys = []
                    for character in line_text:
                        if character == '?':
                            replace_keys.append(index)
                        index += 1
                    replace_count = len(replace_keys) / 2
                    for replacement in range(0, replace_count + 1):
                        string_beginning = replace_keys[replacement]
                        string_end = replace_keys[replacement + 1]
                        replacement_variable_number = int(line_text[string_beginning + 1:string_end])
                        old_string = line_text[string_beginning:string_end + 1]
                        replacement_string = str(design_var[replacement_variable_number + 1])
                        line_text.replace(old_string, replacement_string)
                line_list.append(line_text)
            return line_list

    import shutil
    import os
    files = os.listdir(common_bem_file_folder)
    for bem_input_file in files:
        file = common_bem_file_folder + '/' + bem_input_file
        new_file_text = change_case_file(file, design_vars)
        shutil.copy(file, copy_folder)
        change_file = copy_folder = '/' + bem_input_file
        with open(change_file, 'w') as f:
            for new_string in new_file_text:
                f.write(new_string)
    return


def run_wamit(output_folder, bem_command):
    """
    This function runs the boundary element solver WAMIT in the current case directory.
    :param output_folder:
    :param bem_command:
    :return:
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


input_file_names = parse_input(args.input)
cases_file = input_file_names['cases_file']
design_data = np.genfromtxt(cases_file, delimiter=',')
design_count = design_data.shape[0]

for case in range(design_count):
    print(case)
    design_variables = design_data[case, :]
    print(design_variables)
    case_output_folder = create_case_directory(case + 1, input_file_names['output_directory'])
    create_case_files(input_file_names['common_file_directory'], case_output_folder, design_variables)
    if args.run:
        run_wamit(case_output_folder, input_file_names['run_wamit_directory'])
        read_output()

# if __name__ == "__main__":
#    main()
