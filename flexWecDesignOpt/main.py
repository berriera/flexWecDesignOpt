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


def parse_input():
    """

    :return:
    """
    import yaml
    with open('C:/Users/13365/Documents/GitHub/flexWecDesignOpt/flexWecDesignOpt/examples/barge/input.yaml', 'r') as f:
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


def create_case_files(common_bem_file_folder, copy_folder):
    """

    :param common_bem_file_folder:
    :param copy_folder:
    :return:
    """
    import shutil
    import os
    files = os.listdir(common_bem_file_folder)
    for bem_input_file in files:
        file = common_bem_file_folder + '/' + bem_input_file
        shutil.copy(file, copy_folder)
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
    subprocess.Popen('powershell.exe ' + bem_command)
    # subprocess.run(bem_command)
    return


def read_output():
    """

    :return:
    """
    return


input_file_names = parse_input()
cases_file = input_file_names['cases_file']
design_data = np.genfromtxt(cases_file, delimiter=',')
design_count = design_data.shape[0]

for case in range(design_count):
    print(case)
    design_variables = design_data[case, :]
    print(design_variables)
    case_output_folder = create_case_directory(case + 1, input_file_names['output_directory'])
    create_case_files(input_file_names['common_file_directory'], case_output_folder)
    # run_wamit(case_output_folder, input_file_names['run_wamit_directory'])
    read_output()

# if __name__ == "__main__":
#    main()
