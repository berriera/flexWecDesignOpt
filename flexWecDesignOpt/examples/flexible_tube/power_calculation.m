function P_annual = power_calculation(file_folder, froude_scale, Xs, rs, ro, n, rho, ls, tube_length)

    % Load mode information array generated from Python
    radial_displacement_file_location = fullfile(file_folder, 'radial_displacement.mat');
    load(radial_displacement_file_location);
    f = radial_displacement_array;

    % Hardcoded values used for testing
    %Xs = 0.005; tube_length = 10; rs = 0.274; ls = 1; rho = 1000; n = 0.323; ro = 0.2108;

    % Calculate tube area and array of x values along the tube length
    Ss = pi * (rs ^ 2);
    x = linspace(-tube_length / 2, tube_length / 2, size(f, 2));

    % Read in WAMIT data from output file
    output_file_location = fullfile(file_folder, 'tube.out')
    hydro = struct();
    hydro = Read_WAMIT(hydro, output_file_location,'rao');
    re = squeeze(hydro.ex_re);
    ma = squeeze(hydro.ex_ma);
    ph = squeeze(hydro.ex_ph);
    T = hydro.T;
    w = hydro.w;
    clear hydro;

    % Sum modal responses
    N = length(ma(:,1))-6;
    RAO = zeros(length(T),length(x));  % Based on phasor addition (this is probably correct, based on barge and column results)
    for k = 1:length(T)
        for j = 7:(6+N)
            for i = 7:(6+N)
                RAO(k,:) = RAO(k,:)+...
                    (Xs*ma(i,k)*f(i-6,:)).*(Xs*ma(j,k)*f(j-6,:)).*cos(ph(i,k)-ph(j,k));
            end
        end
    end
    RAO = RAO.^(1/2);
    dL = RAO/rs;

    % Calculate scaled power
    p = load('C:\Users\13365\Desktop\T_Prob_Dist.mat'); % Humbolt Bay Probabilities
    Pa = interp1(p.Ta/sqrt(froude_scale),p.Pa/100,T,'nearest','extrap'); % Assuming exp model is 1:froude_scale scale
    P_annual = 0;
    H = 0.2*ls;   % 0.2 m wave height
    Y = (H/2)*dL; % Strain amplitude
    for k = 1:length(T)
        P(k) = (8*rho*n*(rs/ro)*pi*pi*Ss/T(k)/T(k))*trapz(x,Y(k,:).^2);
        P_annual = P_annual+P(k)*Pa(k);
    end
    P_annual = P_annual*(froude_scale^(7/2))/1000;  % Value at Full Scale (kW), paper approx 100-200 kW
end
