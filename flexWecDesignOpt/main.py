import argparse
import sys
import os
from file_mgmt import parse_input
from file_mgmt import create_case_directory
from file_mgmt import current_case_directory
from substitution import create_case_files
from analysis import run_wamit
from output import read_output
from mesh import create_mesh_file
from mesh import submerged_mesh
from device_types import *

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
                    required=True,
                    help='file path to input.yaml file'
                    )
parser.add_argument('-m', '--mesh',
                    action='store_true',
                    default=False,
                    help='bool: option to use pygmsh and GMSH to parameterize mesh'
                    )
parser.add_argument('-r', '--run',
                    action='store_true',
                    default=False,
                    help='bool: option to run boundary element method command on created input files'
                    )
args = parser.parse_args(sys.argv[1:])


# TODO: add in test mesh refinement argparse option to pop in first design and visualize it in meshmagick,
#   with keyboard option to keep going, retry with smaller mesh refinement factor, or cancel


def main():
    # Grabs information from inputted .yaml file
    input_file_names = parse_input(args.input)
    device_name = input_file_names['device_name']
    common_file_directory = os.path.abspath(input_file_names['common_file_directory'])
    cases_file = os.path.abspath(input_file_names['cases_file'])
    output_directory = os.path.abspath(input_file_names['output_directory'])
    # TODO: create output directory if it doesn't exist
    if args.run:
        try:
            bem_command = input_file_names['run_wamit_command']
        except KeyError:
            print("Boundary element command not included in input file")

    if args.mesh:
        try:
            gmsh_exe_location = os.path.abspath(input_file_names['gmsh_exe_location'])
        except KeyError:
            print("Check find GMSH executable file for meshing.")
        try:
            mesh_refinement_factor = float(input_file_names['mesh_refinement_factor'])
        except KeyError:
            print("Cannot find given mesh refinement factor. Using a default of 0.50.")

    # Creates array of all design variables from inputted .csv file
    design_data = np.genfromtxt(cases_file, delimiter=',')
    design_count = design_data.shape[0]

    # Finds corresponding class name according to the inputted device name
    device_type = getattr(sys.modules[__name__], device_name)

    for case in range(design_count):
        print('Case:', str(case + 1))
        design_variables = design_data[case, :]

        # Obtains device class information
        device = device_type(design_variables)
        substitution_array = device.substitutions()
        geometry = device.geometry()

        # Creates case directory then copies and changes boundary element method input files into the directory
        case_output_folder = create_case_directory(output_directory, case + 1)
        # case_output_folder = current_case_directory()
        create_case_files(common_file_directory, case_output_folder, substitution_array)

        if args.mesh:
            print('\tMeshing...')
            create_mesh_file(geometry, device_name, case_output_folder,
                                           gmsh_exe_location, mesh_refinement_factor)
            submerged_mesh(device_name)

        if args.run:
            print('\tRunning BEM...')
            try:
                run_wamit(case_output_folder, bem_command)
            except UnboundLocalError:
                print('\tCannot find boundary element model command in input file.')
            except FileNotFoundError:
                print('\tCannot run given boundary element model command.')
            read_output()
        print('\tDone.')

    print('\n\nAll cases completed.')


if __name__ == "__main__":
    main()
