% Four-Tank System: Step Response Simulation and Identification

current_dir = fileparts(mfilename('fullpath'));

%   Add common_functions folder to PATH
folder_to_add = fullfile(current_dir, 'common_functions');
addpath(folder_to_add);
disp(['Added folder to path: ', folder_to_add]);

% --------------------------------------------------------------
% Parameters (from original script)
% --------------------------------------------------------------

a1 = 1.2272;        %[cm2] Area of outlet pipe 1
a2 = 1.2272;        %[cm2] Area of outlet pipe 2
a3 = 1.2272;        %[cm2] Area of outlet pipe 3
a4 = 1.2272;        %[cm2] Area of outlet pipe 4
A1 = 380.1327;      %[cm2] Cross sectional area of tank 1
A2 = 380.1327;      %[cm2] Cross sectional area of tank 2
A3 = 380.1327;      %[cm2] Cross sectional area of tank 3
A4 = 380.1327;      %[cm2] Cross sectional area of tank 4
gamma1 = 0.45;      % Flow distribution constant. Valve 1
gamma2 = 0.40;      % Flow distribution constant. Valve 2
g = 981;            % [cm/s2] The acceleration of gravity
rho = 1.00;         % [g/cm3] Density of water
p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];
A = [A1; A2; A3; A4];

% --------------------------------------------------------------

% Simulation scenario

t0 = 0.0;                       % [s] Initial time
t_final = 2*60;                 % [s] Final time (2 min)
m10 = 0.0;                      % [g] Liquid mass in tank 1 at time t0
m20 = 0.0;                      % [g] Liquid mass in tank 2 at time t0
m30 = 0.0;                      % [g] Liquid mass in tank 3 at time t0
m40 = 0.0;                      % [g] Liquid mass in tank 4 at time t0
F1 = 0;                         % [cm3/s] Flow rate from pump 1
F2 = 0;                         % [cm3/s] Flow rate from pump 2
x0 = [m10; m20; m30; m40];
u = [F1; F2];                   % Initially set to 0

% --------------------------------------------------------------

dt = 0.001;
t = t0:dt:t_final;
N = length(t);
norm_d = 50;
Rd = [0.2 0.1 0.1 0.2];

%   Scalar wiener noise generation and iteration index
[~,~,dw] = ScalarStdWienerProcess(t(end),N,2,5);
dw_it = 1;

d_diff = [0.0; 0.0];
d_mean = [norm_d/2;norm_d/2];           %   Set to the mean value of d
