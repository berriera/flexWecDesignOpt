
    output_file_location = fullfile(file_folder, 'tube.out')
    hydro = struct();
    hydro = Read_WAMIT(hydro, output_file_location,'rao');
    re = squeeze(hydro.ex_re);
    ma = squeeze(hydro.ex_ma);
    ph = squeeze(hydro.ex_ph);
    T = hydro.T;
    w = hydro.w;