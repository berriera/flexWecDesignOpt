import argparse
import sys
import numpy as np
import string
from parse_input import parse_input
from create_case_directory import create_case_directory
from create_case_files import create_case_files
from run_wamit import run_wamit
from read_output import read_output
from create_mesh_file_from_geometry import create_mesh_file_from_geometry
from device_types import *

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
                    required=True,
                    help='file path to input.yaml file'
                    )
parser.add_argument('-r', '--run',
                    action='store_true',
                    default=False,
                    help='bool: option to run boundary element method command on created input files'
                    )
args = parser.parse_args(sys.argv[1:])


# TODO: add in test mesh refinement argparse option to pop in first design and visualize it in meshmagick,
#   with keyboard option to keep going, retry with smaller mesh refinement factor, or cancel

# TODO: add in option to pass in array; if this happens, read it in as a numpy array

def main():
    # Grabs information from inputted .yaml file
    input_file_names = parse_input(args.input)
    device_name = input_file_names['device_name']
    common_file_directory = input_file_names['common_file_directory']
    cases_file = input_file_names['cases_file']
    output_directory = input_file_names['output_directory']
    bem_command = input_file_names['run_wamit_directory']
    gmsh_exe_location = input_file_names['gmsh_exe_location']

    # Creates array of all design variables from inputted .csv file
    design_data = np.genfromtxt(cases_file, delimiter=',')
    design_count = design_data.shape[0]

    # Finds corresponding class name according to the inputted device name
    device_type = getattr(sys.modules[__name__], string.capwords(device_name))

    for case in range(design_count):
        print('Case:', str(case + 1))
        design_variables = design_data[case, :]

        # Obtains device class information
        device = device_type(design_variables)
        substitution_array = device.substitutions()
        geometry = device.geometry()

        # Creates case directory then copies and changes boundary element method input files into the directory
        case_output_folder = create_case_directory(case + 1, output_directory)
        create_case_files(common_file_directory, case_output_folder, substitution_array)

        print('\tMeshing...')
        create_mesh_file_from_geometry(geometry, device_name, case_output_folder, gmsh_exe_location)

        if args.run:
            print('\tRunning BEM...')
            run_wamit(case_output_folder, bem_command)
            read_output()
        print('\tDone.')

    print('\n\nAll cases completed.')


if __name__ == "__main__":
    main()
