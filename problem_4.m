clc; clear; close all
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
tf = 60*30;      % [s]
Ts = 1;          % [s]
t = t0:Ts:tf;
N = length(t);

F1 = 300; F2 = 300;
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

%% Step responses fot deterministic model
% --- Define step changes in references ---

outputFolder = false;
u = [F1*1.1;F2]
%determinitic simulation step 10% in F1
% outputFolder = fullfile('figures', 'problem_4', 'step_10_F1','deterministic'); 
% if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
[t, y_10_F1] = simulate_deterministic(t, xs, u, d, p, false, outputFolder);


u = [F1*1.25; F2]
%determinitic simulation step 25% in F1
% outputFolder = fullfile('figures', 'problem_4', 'step_25_F1','deterministic'); 
% if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
[t, y_25_F1] = simulate_deterministic(t, xs, u, d, p, false, outputFolder);

u = [F1*1.5; F2]
%determinitic simulation step 25% in F1
% outputFolder = fullfile('figures', 'problem_4', 'step_50_F1','deterministic'); 
% if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
[t, y_50_F1] = simulate_deterministic(t, xs, u, d, p, false, outputFolder);

% Assuming y_10_F1, y_25_F1, y_50_F1 are 4×N matrices of heights (or masses)
[fig, ax] = plot_step_responses(t, y_10_F1, y_25_F1, y_50_F1, p, "F_1");
outputFolder = fullfile('figures', 'problem_4', 'steps_F1','deterministic'); 

if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(fig, fullfile(outputFolder, 'Steps_F1.pdf'), 'ContentType', 'vector');

outputFolder = false;
u = [F1;F2*1.1]
%determinitic simulation step 10% in F2
% outputFolder = fullfile('figures', 'problem_4', 'step_10_F1','deterministic'); 
% if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
[t, y_10_F2] = simulate_deterministic(t, xs, u, d, p, false, outputFolder);

u = [F1; F2*1.25]
%determinitic simulation step 25% in F2
% outputFolder = fullfile('figures', 'problem_4', 'step_25_F1','deterministic'); 
% if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
[t, y_25_F2] = simulate_deterministic(t, xs, u, d, p, false, outputFolder);

u = [F1; F2*1.5]
%determinitic simulation step 25% in F2
% outputFolder = fullfile('figures', 'problem_4', 'step_50_F1','deterministic'); 
% if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
[t, y_50_F2] = simulate_deterministic(t, xs, u, d, p, false, outputFolder);

% Assuming y_10_F1, y_25_F1, y_50_F1 are 4×N matrices of heights (or masses)
[fig, ax] = plot_step_responses(t, y_10_F2, y_25_F2, y_50_F2, p, "F_2");
outputFolder = fullfile('figures', 'problem_4', 'steps_F2','deterministic'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(fig, fullfile(outputFolder, 'Steps_F2.pdf'), 'ContentType', 'vector');


%% Normalized step responses for F1
% y_10_F1, y_25_F1, y_50_F1  → 4×N matrices
% hs → steady-state heights [4x1]
% F1, F2, u_s already defined

% Define input steps
du_10_F1 = (F1*1.1 - F1);    % step magnitude in cm³/s
du_25_F1 = (F1*1.25 - F1);
du_50_F1 = (F1*1.5 - F1);

% Normalize outputs
y_norm_10_F1 = (y_10_F1 - hs) / du_10_F1;
y_norm_25_F1 = (y_25_F1 - hs) / du_25_F1;
y_norm_50_F1 = (y_50_F1 - hs) / du_50_F1;

% Plot normalized step responses
[fig_norm_F1, ax_norm_F1] = plot_step_responses(t, y_norm_10_F1, y_norm_25_F1, y_norm_50_F1, p, 'F_1');
sgtitle('Normalized Step Responses to F_1 Steps', 'FontSize', 14);
outputFolder = fullfile('figures', 'problem_4', 'normalized_steps','deterministic'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(fig_norm_F1, fullfile(outputFolder, 'normalized_steps_F1.pdf'), 'ContentType', 'vector');

%Normalized steps for F2
du_10_F2 = (F2*1.1 - F2);
du_25_F2 = (F2*1.25 - F2);
du_50_F2 = (F2*1.5 - F2);

y_norm_10_F2 = (y_10_F2 - hs) / du_10_F2;
y_norm_25_F2 = (y_25_F2 - hs) / du_25_F2;
y_norm_50_F2 = (y_50_F2 - hs) / du_50_F2;

[fig_norm_F2, ax_norm_F2] = plot_step_responses(t, y_norm_10_F2, y_norm_25_F2, y_norm_50_F2, p, 'F_2');
sgtitle('Normalized Step Responses to F_2 Steps', 'FontSize', 14);
exportgraphics(fig_norm_F2, fullfile(outputFolder, 'normalized_steps_F2.pdf'), 'ContentType', 'vector');

du_10_F1 = F1*1.1 - F1;
y_norm_F1 = (y_10_F1 - hs) / du_10_F1;

du_10_F2 = F2*1.1 - F2;
y_norm_F2 = (y_10_F2 - hs) / du_10_F2;

[fig_tf, ax_tf] = plot_transfer_identification(t, y_norm_F1, y_norm_F2);
outputFolder = fullfile('figures', 'problem_4', 'normalized_tf','deterministic'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(fig_tf, fullfile(outputFolder, 'normalized_tf.pdf'), 'ContentType', 'vector');

% % 
%%
disp('=== IDENTIFYING TRANSFER FUNCTIONS WITH TFEST ===');

Ts = 1; % [s] sampling time

% --- Input signal: unit step (same length as t)
u_step = ones(length(t),1);

% --- Create iddata objects for each input/output pair
data_G11 = iddata(y_norm_F1(1,:)', u_step, Ts);
data_G21 = iddata(y_norm_F1(2,:)', u_step, Ts);
data_G12 = iddata(y_norm_F2(1,:)', u_step, Ts);
data_G22 = iddata(y_norm_F2(2,:)', u_step, Ts);

% --- Model orders
np_direct = 1;   % 1 pole for G11, G22
np_cross  = 2;   % 2 poles for G12, G21
nz = 0;          % 0 zeros (strictly proper)

% --- Identification
G11_est = tfest(data_G11, np_direct, nz);
G22_est = tfest(data_G22, np_direct, nz);
G12_est = tfest(data_G12, np_cross,  nz);
G21_est = tfest(data_G21, np_cross,  nz);

% --- Display results
disp('Identified transfer functions (tfest results):');
G11_est
G12_est
G21_est
G22_est

% --- Compare fits
figure('Name','TFEST Model Fits');
subplot(2,2,1); compare(data_G11, G11_est); title('G_{11}: h_1 \leftarrow F_1');
subplot(2,2,2); compare(data_G12, G12_est); title('G_{12}: h_1 \leftarrow F_2');
subplot(2,2,3); compare(data_G21, G21_est); title('G_{21}: h_2 \leftarrow F_1');
subplot(2,2,4); compare(data_G22, G22_est); title('G_{22}: h_2 \leftarrow F_2');
sgtitle('Transfer Function Identification using TFEST');

%%
fprintf('\n--- Extracted Parameters from TFEST ---\n');

% Helper function to extract parameters from any tf
getParams = @(sys) deal(dcgain(sys), -1./real(pole(sys)));

[K11_est, tau11] = getParams(G11_est);
[K22_est, tau22] = getParams(G22_est);
[K12_est, tau12] = getParams(G12_est);
[K21_est, tau21] = getParams(G21_est);

fprintf('G11(s): K = %.4e, tau = %.2f s\n', K11_est, tau11);
fprintf('G22(s): K = %.4e, tau = %.2f s\n', K22_est, tau22);
fprintf('G12(s): K = %.4e, taus = [%.2f, %.2f] s\n', K12_est, sort(tau12));
fprintf('G21(s): K = %.4e, taus = [%.2f, %.2f] s\n', K21_est, sort(tau21));


% ------------------------------------------------------------
% Supporting residual function
% ------------------------------------------------------------
function r = FourTankSystemSteadyResidual(x, u, d, p)
r = FourTankSystemModified(0, x, u, d, p); % residual f(x)=0
end
