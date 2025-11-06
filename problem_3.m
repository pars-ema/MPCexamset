%% PID CONTROL of the SDE Nonlinear Four-Tank System
clear; close all; clc;
% 
% set(groot, 'defaultTextInterpreter', 'latex');
% set(groot, 'defaultLegendInterpreter', 'latex');
% set(groot, 'defaultAxesTickLabelInterpreter', 'latex');


% -----------------------------------------------------------
% Parameters
% -----------------------------------------------------------
a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;        % [cm2]
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327; % [cm2]
g = 981; rho = 1.00;
gamma1 = 0.58; gamma2 = 0.72;
p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2; rho];

% ------------------------------------------------------------
% Simulation scenario
% ------------------------------------------------------------
t0 = 0.0;        % [s]
tf = 60*60;      % [s]
Ts = 1;          % [s]
t = t0:Ts:tf;
N = length(t);

% Mean disturbances
F3 = 250; F4 = 250;           % [cm3/s]
d = [F3; F4];

% Steady-state inputs (nominal pump flows)
u_s = [300; 300];             % [cm3/s]

% ------------------------------------------------------------
% Compute steady state automatically using fsolve
% ------------------------------------------------------------
xs0 = 500*ones(4,1);  % initial guess
opts = optimoptions('fsolve','Display','none','FunctionTolerance',1e-12,'StepTolerance',1e-12);
xs = fsolve(@(x) FourTankSystemSteadyResidual(x,u_s,d,p), xs0, opts);

% Compute steady-state heights for reference
hs = xs ./ (rho * [A1; A2; A3; A4]);

disp('Steady-state heights [cm]:'); disp(hs.');
disp('Steady-state masses [g]:'); disp(xs.');

%%

close all

% --- Steady-state heights (from fsolve) ---
hs

% --- Define step changes in references ---
t_change_h1 = [0 600 2400];         % [s]
values_h1   = [hs(1) hs(1)*1.2 hs(1)*0.8];

t_change_h2 = [0 900 3000];         % [s]
values_h2   = [hs(2) hs(2) hs(2)];

% --- Generate reference trajectories ---
[r1, r2] = generate_stair_references(t, t_change_h1, values_h1, t_change_h2, values_h2);
h_ref_total = [r1; r2];

% Plot the references
figure;
stairs(t, r1, 'LineWidth', 1.5, 'DisplayName', 'Tank 1 Reference');
hold on;
stairs(t, r2, 'LineWidth', 1.5, 'DisplayName', 'Tank 2 Reference');
xlabel('Time [s]');
ylabel('Reference Level [cm]');
title('Piecewise Constant Tank Level References');
legend('Location','best');
grid on;


%% Only P control fo deterministic
outputFolder = fullfile('figures', 'problem_3', 'P', 'deterministic');
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
% Call the function
Kp = [2.5;2.5];
Ki = [0;0];
Kd = [0;0];
simulate_pid_deterministic(t, xs, u_s, d, p, h_ref_total, Kp, Ki, Kd, outputFolder);
%%
% srochastic
outputFolder = fullfile('figures', 'problem_3', 'P', 'stochastic');
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
% Call the function
simulate_pid_stochastic(t, xs, u_s, p, h_ref_total,Kp, Ki, Kd, outputFolder);

%% PI control
outputFolder = fullfile('figures', 'problem_3', 'PI', 'deterministic');
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
% Call the function
Kp = [2.5;2.5];
Ki = [0.1;0.1];
Kd = [0;0];
simulate_pid_deterministic(t, xs, u_s, d, p, h_ref_total, Kp, Ki, Kd, outputFolder);


% stochastic
outputFolder = fullfile('figures', 'problem_3', 'PI', 'stochastic');
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
% Call the function
simulate_pid_stochastic(t, xs, u_s, p, h_ref_total,Kp, Ki, Kd, outputFolder);

%% PID
outputFolder = fullfile('figures', 'problem_3', 'PID', 'deterministic');
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
% Call the function
Kp = [2.5;2.5];
Ki = [0.1;0.1];
Kd = [0.1;0.1];
simulate_pid_deterministic(t, xs, u_s, d, p, h_ref_total, Kp, Ki, Kd, outputFolder);


% stochastic
outputFolder = fullfile('figures', 'problem_3', 'PID', 'stochastic');
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
% Call the function
simulate_pid_stochastic(t, xs, u_s, p, h_ref_total,Kp, Ki, Kd, outputFolder);


%% More compex experiment 
% --- Define step changes in references ---
t_change_h1 = [0 600 2400];         % [s]
values_h1   = [hs(1) hs(1)*1.2 hs(1)*0.8];

t_change_h2 = [0 900 3000];         % [s]
values_h2   = [hs(2) hs(2)*0.7 hs(2)*1.3];

% --- Generate reference trajectories ---
[r1, r2] = generate_stair_references(t, t_change_h1, values_h1, t_change_h2, values_h2);
h_ref_total = [r1; r2];

% Plot the references
figure;
stairs(t, r1, 'LineWidth', 1.5, 'DisplayName', 'Tank 1 Reference');
hold on;
stairs(t, r2, 'LineWidth', 1.5, 'DisplayName', 'Tank 2 Reference');
xlabel('Time [s]');
ylabel('Reference Level [cm]');
title('Piecewise Constant Tank Level References');
legend('Location','best');
grid on;

%PID
outputFolder = fullfile('figures', 'problem_3', 'PID_2', 'deterministic');
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
% Call the function
Kp = [2.5;2.5];
Ki = [0.1;0.1];
Kd = [0.1;0.1];
simulate_pid_deterministic(t, xs, u_s, d, p, h_ref_total, Kp, Ki, Kd, outputFolder);


% stochastic
outputFolder = fullfile('figures', 'problem_3', 'PID_2', 'stochastic');
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
% Call the function
simulate_pid_stochastic(t, xs, u_s, p, h_ref_total,Kp, Ki, Kd, outputFolder);



% ------------------------------------------------------------
% Supporting residual function
% ------------------------------------------------------------
function r = FourTankSystemSteadyResidual(x, u, d, p)
r = FourTankSystemModified(0, x, u, d, p); % residual f(x)=0
end
