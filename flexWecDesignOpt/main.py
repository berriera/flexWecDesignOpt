import argparse
import sys
import numpy as np
from parse_input import parse_input
from create_case_directory import create_case_directory
from create_case_files import create_case_files
from run_wamit import run_wamit
from read_output import read_output

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
