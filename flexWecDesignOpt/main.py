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
    import yaml
    with open('C:/Users/13365/Documents/GitHub/flexWecDesignOpt/flexWecDesignOpt/examples/barge/input.yaml', 'r') as f:
        file_names = yaml.safe_load(f)
    return file_names


def create_case_directory(case_number, output_directory):
    import os
    path = output_directory + '/case_' + str(case_number)
    os.makedirs(path)


def create_case_files():
    return


def run_wamit():
    # subprocess.run(["C:\WAMITv7", "runwamit", "test01"]) # can also try subprocess.Popen
    return


def read_output():
    return


input_file_names = parse_input()
cases_file = input_file_names['cases_file']
design_data = np.genfromtxt(cases_file, delimiter=',')
design_count = design_data.shape[0]

for case in range(1, design_count + 1):
    print(case)
    design_variables = design_data[case, :]
    print(design_variables)
    create_case_directory(case, input_file_names['output_directory'])
    create_case_files()
    run_wamit()
    read_output()

# if __name__ == "__main__":
#    main()
