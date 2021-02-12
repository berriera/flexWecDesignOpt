function read_WAMIT_output(output_file_name)

    % Change directory to location of output file
    [case_directory, ~, ~] = fileparts(output_file_name)
    cd(case_directory)

    % Fill results struct and save to .mat for processing to dict in Python
    hydro = struct();
    hydro = Read_WAMIT(hydro, output_file_name,'rao');
    clear output_file_name case_directory
    save('WAMIT_results.mat');
end
