import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True,
                    help='file path to input.yaml file'
                    )
parser.add_argument('-r', '--run',
                    default=True,
                    type=bool,
                    help='run boundary element method software on created input files'
                    )

args = parser.parse_args(sys.argv[1:])


case_count = parse_input()

for case in case_count:
    case_data = parse_input(case)
    create_case_directory()
    create_case_files()
    run_wamit()
    read_output()

if __name__ == "__main__":
    main()