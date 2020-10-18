clc; clear all; close all;

% Load variables from Python
load('run.mat')
load('radial_displacement_function.mat')

% Maximum modal radial displacement variable
Xs = 0.005;

% Read in WAMIT data from output file
hydro = struct();
hydro = Read_WAMIT(hydro,'tube.out','rao');
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

% Power Calculation
p = load('T_Prob_Dist.mat'); % Humbolt Bay Probabilities
Pa = interp1(p.Ta/sqrt(10),p.Pa/100,T,'nearest','extrap'); % Assuming exp model is 1:10 scale
Pannual = 0;
H = 0.2*ls;
Y = (H/2)*dL; % Strain amplitude
for k = 1:length(T)
    P(k) = (8*rho*n*(rs/ro)*pi*pi*Ss/T(k)/T(k))*trapz(x,Y(k,:).^2);
    Pannual = Pannual+P(k)*Pa(k);
end
Pannual = Pannual*(10^(7/2))/1000  % Value at Full Scale (kW), paper approx 100-200 kW