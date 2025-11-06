clear; close all; clc;

% -----------------------------------------------------------
% Parameters
% -----------------------------------------------------------
a1 = 1.2272 %[cm2] Area of outlet pipe 1
a2 = 1.2272 %[cm2] Area of outlet pipe 2
a3 = 1.2272 %[cm2] Area of outlet pipe 3
a4 = 1.2272 %[cm2] Area of outlet pipe 4
A1 = 380.1327 %[cm2] Cross sectional area of tank 1
A2 = 380.1327 %[cm2] Cross sectional area of tank 2
A3 = 380.1327 %[cm2] Cross sectional area of tank 3
A4 = 380.1327 %[cm2] Cross sectional area of tank 4
g = 981 %[cm/s2] The acceleration of gravity
rho = 1.00;
gamma1 = 0.58; % Flow distribution constant. Valve 1
gamma2 = 0.72; % Flow distribution constant. Valve 2
p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2; rho]


% ------------------------------------------------------------
% Simulation scenario
% ------------------------------------------------------------
t0 = 0.0; % [s] Initial time
tf = 60*60; % [s] Final time
m10 = 0.0; % [g] Liquid mass in tank 1 at time t0
m20 = 0.0; % [g] Liquid mass in tank 2 at time t0
m30 = 0.0; % [g] Liquid mass in tank 3 at time t0
m40 = 0.0; % [g] Liquid mass in tank 4 at time t0
F1 = 300; % [cm3/s] Flow rate from pump 1
F2 = 300; % [cm3/s] Flow rate from pump 2
x0 = [m10; m20; m30; m40];
u = [F1; F2];
F3 = 250;
F4 = 250;
d = [F3; F4];

% Time setup
t0 = 0;
tf = 60*60;
Ts = 1;
t  = t0:Ts:tf;


%% problem 2.1 Determinitic nonlinear model with constant input u
%determinitic simulation with constant input and constant disturbances
outputFolder = fullfile('figures', 'problem_2', 'deterministic'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
simulate_deterministic(t, x0, u, d, p, true, outputFolder);

%% problem 2.2 Stochastic Nonlinear model
%stochastic nonlinear simulation with constant input and piecewise constant disturbances
outputFolder = fullfile('figures', 'problem_2', 'stochastic'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
simulate_stochastic(t, x0, u, p, outputFolder);

%% problrm 2.3 SDE Non-linear model
%sde simulation
outputFolder = fullfile('figures', 'problem_2', 'sde'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
simulate_sde(t, x0, u, d, p, outputFolder);


%% Problem 2.4 simulation with random step inputs
[u1, u2, tc1, tc2, v1, v2] = generate_random_stairs_inputs(t, 5, 0, 500, 42);
u = [u1; u2];  % Combine into 2xN input matrix

% Disturbances
d = [250; 250];
% Initial condition
x0 = zeros(4,1);


%determinitic simulation with constant input and constant disturbances
outputFolder = fullfile('figures', 'problem_2', 'simulation','deterministic'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
simulate_deterministic(t, x0, u, d, p, true, outputFolder);

%stochastic nonlinear simulation with constant input and piecewise constant disturbances
outputFolder = fullfile('figures', 'problem_2', 'simulation', 'stochastic'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
simulate_stochastic(t, x0, u, p, outputFolder);

%sde simulation
outputFolder = fullfile('figures', 'problem_2', 'simulation', 'sde'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
simulate_sde(t, x0, u, d, p, outputFolder);


