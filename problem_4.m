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
t_final = 60*30;      % [s]
Ts = 1;          % [s]
t = t0:Ts:t_final;
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

% --- Add annotations AFTER figure creation ---
% Access the desired axes (for example, G11 in position (1,1))
axes(ax_tf(1,1));  % or use 'nexttile' if using tiledlayout
% Add lines and text
xline(198, 'r--', 'LineWidth', 1.2);
text(198, 0.5, '$T_1$', 'Interpreter','latex', ...
     'Rotation',90, 'VerticalAlignment','bottom', ...
     'HorizontalAlignment','right', 'Color','red');

yline(0.632*0.301, 'r--', 'LineWidth', 1.2);
text(1400, 0.632*0.301, '$0.632 \cdot K_{11}$', 'Interpreter','latex', ...
     'VerticalAlignment','bottom', 'Color','red');

axes(ax_tf(2,2));  % G22 plot
xline(220, 'r--', 'LineWidth', 1.2);
text(220, 0.5, '$T_2$', 'Interpreter','latex', ...
     'Rotation',90, 'VerticalAlignment','bottom', ...
     'HorizontalAlignment','right', 'Color','red');
yline(0.632*0.415, 'r--', 'LineWidth', 1.2);
text(1400, 0.632*0.415, '$0.632\cdot K_{22}$', 'Interpreter','latex', ...
     'VerticalAlignment','bottom', 'Color','red');
exportgraphics(fig_tf, fullfile(outputFolder, 'normalized_tf_with_identification.pdf'), 'ContentType', 'vector');




%% Normalized step responses for determinitic model with 10% step
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
outputFolder = fullfile('figures', 'problem_4', 'normalized_steps','sde'); 
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
outputFolder = fullfile('figures', 'problem_4', 'normalized_tf','sde'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(fig_tf, fullfile(outputFolder, 'normalized_tf.pdf'), 'ContentType', 'vector');


%% 4.3 Identifying transfer function parameters
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

%%
disp('--- Computing discrete-time Markov parameters using mimoctf2dss ---');

% --- Sampling settings ---
Ts   = 10;      % [s] sampling time
Nmax = 100;    % maximum state dimension
tol  = 1e-8;   % numerical tolerance

% --- Extract numerator and denominator directly from identified models ---
num = cell(2,2);
den = cell(2,2);
lambda = zeros(2,2);

num{1,1} = G11_est.Numerator;  den{1,1} = G11_est.Denominator;
num{1,2} = G12_est.Numerator;  den{1,2} = G12_est.Denominator;
num{2,1} = G21_est.Numerator;  den{2,1} = G21_est.Denominator;
num{2,2} = G22_est.Numerator;  den{2,2} = G22_est.Denominator;

% --- Convert continuous-time MIMO TF to discrete-time state-space ---
[Ad, Bd, Cd, Dd, sH] = mimoctf2dss(num, den, lambda, Ts, Nmax, tol);

% Save them to a .mat file
save('discrete_model_from_step.mat', 'Ad', 'Bd', 'Cd', 'Dd', 'sH');


% --- Compute Markov parameters ---
Nimp = 1800/Ts;
ny = size(Cd,1);
nu = size(Bd,2);
H = zeros(ny, nu, Nimp);
H(:,:,1) = Dd;
for k = 2:Nimp
    H(:,:,k) = Cd * (Ad^(k-2)) * Bd;
end

% --- Save the Markov parameters for later comparison ---
H_identified_problem_4 = H;  % save in a clear variable
save('H_identified_problem_4.mat','H_identified_problem_4','Ts');
disp('✅ Saved identified Markov parameters to Markov_identified.mat');


% --- Plot Markov parameters ---    
figure('Name','Discrete-Time Markov Parameters');
titles = {'$h_{11}[k]$', '$h_{12}[k]$', '$h_{21}[k]$', '$h_{22}[k]$'};
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j);
        plot(0:Nimp-1, squeeze(H(i,j,:)), 'rx',  'LineWidth', 1.3);
        grid on;
        xlabel('Sample $k$', 'Interpreter','latex', 'FontSize',13);
        ylabel('Amplitude', 'FontSize',12);
        title(titles{(i-1)*2+j}, 'Interpreter','latex', 'FontSize',14);
    end
end
sgtitle('Impulse Response Coefficients (Markov Parameters)', ...
        'Interpreter','latex','FontSize',15);

% --- Optional: Hankel singular values ---
figure('Name','Hankel Singular Values');
semilogy(sH,'b.-','LineWidth',1.5,'MarkerSize',10);
grid on;
xlabel('State number');
ylabel('Singular Value');
title('Hankel Singular Values (Model Order Content)', 'FontSize',12);

disp('✅ Markov parameter computation complete.');


%% suggestion from chat
%% === Impulse responses of estimated transfer functions ===
disp('--- Plotting impulse responses for estimated transfer functions ---');

% Convert identified idtf models to tf
G11_tf = tf(G11_est);
G12_tf = tf(G12_est);
G21_tf = tf(G21_est);
G22_tf = tf(G22_est);

titles = {
    '$G_{11}$ : $h_1 \leftarrow F_1$', ...
    '$G_{12}$ : $h_1 \leftarrow F_2$', ...
    '$G_{21}$ : $h_2 \leftarrow F_1$', ...
    '$G_{22}$ : $h_2 \leftarrow F_2$'};

% Colors
color_ct = [0 0.45 0.74];   % blue for continuous impulse
color_d  = [0.95 0. 0.1]; % red for discrete Markov

figure('Name','Impulse Responses of Estimated Transfer Functions');
set(gcf, 'Position', [100 100 1000 640]);  % [left bottom width height]

% G11
subplot(2,2,1); hold on; grid on;
[Y, T] = impulse(G11_tf, t);
plot(T, Y, 'Color', color_ct, 'LineWidth', 2, 'DisplayName', 'Continuous-time impulse');
plot((0:Nimp-1)*Ts, squeeze(H(1,1,:))./Ts, 'x', ...
     'Color', color_d, 'LineWidth', 1.3, 'DisplayName', 'Discrete Markov');
xlabel('Time [s]', 'FontSize', 12, 'FontWeight','bold');
ylabel('Amplitude', 'FontSize', 12, 'FontWeight','bold');
title(titles{1}, 'Interpreter','latex', 'FontSize',14);
legend('Location','northeast', 'Interpreter','latex');
set(gca, 'FontSize', 11);

% G12
subplot(2,2,2); hold on; grid on;
[Y, T] = impulse(G12_tf, t);
plot(T, Y, 'Color', color_ct, 'LineWidth', 2, 'DisplayName', 'Continuous-time impulse');
plot((0:Nimp-1)*Ts, squeeze(H(1,2,:))./Ts, 'x', ...
     'Color', color_d, 'LineWidth', 1.3, 'DisplayName', 'Discrete Markov');
xlabel('Time [s]', 'FontSize', 12, 'FontWeight','bold');
ylabel('Amplitude', 'FontSize', 12, 'FontWeight','bold');
title(titles{2}, 'Interpreter','latex', 'FontSize',14);
legend('Location','northeast', 'Interpreter','latex');
set(gca, 'FontSize', 11);

% G21
subplot(2,2,3); hold on; grid on;
[Y, T] = impulse(G21_tf, t);
plot(T, Y, 'Color', color_ct, 'LineWidth', 2, 'DisplayName', 'Continuous-time impulse');
plot((0:Nimp-1)*Ts, squeeze(H(2,1,:))./Ts, 'x', ...
     'Color', color_d, 'LineWidth', 1.3, 'DisplayName', 'Discrete Markov');
xlabel('Time [s]', 'FontSize', 12, 'FontWeight','bold');
ylabel('Amplitude', 'FontSize', 12, 'FontWeight','bold');
title(titles{3}, 'Interpreter','latex', 'FontSize',14);
legend('Location','northeast', 'Interpreter','latex');
set(gca, 'FontSize', 11);

% G22
subplot(2,2,4); hold on; grid on;
[Y, T] = impulse(G22_tf, t);
plot(T, Y, 'Color', color_ct, 'LineWidth', 2, 'DisplayName', 'Continuous-time impulse');
plot((0:Nimp-1)*Ts, squeeze(H(2,2,:))./Ts, 'x', ...
     'Color', color_d, 'LineWidth', 1.3, 'DisplayName', 'Discrete Markov');
xlabel('Time [s]', 'FontSize', 12, 'FontWeight','bold');
ylabel('Amplitude', 'FontSize', 12, 'FontWeight','bold');
title(titles{4}, 'Interpreter','latex', 'FontSize',14);
legend('Location','northeast', 'Interpreter','latex');
set(gca, 'FontSize', 11);

% Save figure
outputFolder = fullfile('figures','problem_4','Markov_parameters');
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
exportgraphics(gcf, fullfile(outputFolder,'impulse_responses_comparison.pdf'),'ContentType','vector');

disp('✅ Impulse responses plotted, labeled, and saved.');



%% 4.2 step responses with process and measurement noise
% --- Define measurement noise levels ---
sigma_values = {[1;1;1;1], [2;2;2;2], [5;5;5;5]};
sigma_labels = {'sigma1','sigma2','sigma5'};  % folder suffixes

% --- Loop through each sigma_measurement configuration ---
for k = 1:length(sigma_values)

    sigma_measurement = sigma_values{k};
    label = sigma_labels{k};

    fprintf('\n=== Running SDE step responses for %s ===\n', label);

    %% === Step in F1 ===
    outputFolder = false;

    u = [F1*1.1; F2];
    [t, ~, y_10_F1] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);

    u = [F1*1.25; F2];
    [t, ~, y_25_F1] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);

    u = [F1*1.5; F2];
    [t, ~, y_50_F1] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);

    % Plot results for F1 steps
    [fig_F1, ax_F1] = plot_step_responses(t, y_10_F1, y_25_F1, y_50_F1, p, sprintf('$F_1, (\\sigma_m^2 = %d)$', sigma_measurement(1)));

    % Save figure
    outputFolder = fullfile('figures', 'problem_4', 'steps_F1', ['sde_' label]);
    if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end
    exportgraphics(fig_F1, fullfile(outputFolder, 'Steps_F1.pdf'), 'ContentType', 'vector');

    %% === Step in F2 ===
    outputFolder = false;

    u = [F1; F2*1.1];
    [t, ~, y_10_F2] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);

    u = [F1; F2*1.25];
    [t, ~, y_25_F2] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);

    u = [F1; F2*1.5];
    [t, ~, y_50_F2] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);

    % Plot results for F2 steps
    [fig_F2, ax_F2] = plot_step_responses(t, y_10_F2, y_25_F2, y_50_F2, p, sprintf('$F_2, (\\sigma_m^2 = %d)$', sigma_measurement(1)));

    % Save figure
    outputFolder = fullfile('figures', 'problem_4', 'steps_F2', ['sde_' label]);
    if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end
    exportgraphics(fig_F2, fullfile(outputFolder, 'Steps_F2.pdf'), 'ContentType', 'vector');

end

fprintf('\nAll SDE step-response simulations completed successfully!\n');

%% Normalized step responses for stochastic model model with medium sigma
% y_10_F1, y_25_F1, y_50_F1  → 4×N matrices
% hs → steady-state heights [4x1]
% F1, F2, u_s already defined

outputFolder = false;
sigma_measurement = sigma_values{2}
u = [F1*1.1; F2];
[t, ~, y_10_F1] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);
u = [F1*1.25; F2];
[t, ~, y_25_F1] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);
u = [F1*1.5; F2];
[t, ~, y_50_F1] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);


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
outputFolder = fullfile('figures', 'problem_4', 'normalized_steps','sde'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(fig_norm_F1, fullfile(outputFolder, 'normalized_steps_F1.pdf'), 'ContentType', 'vector');



%Normalized steps for F2
u = [F1; F2*1.1];
[t, ~, y_10_F2] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);
u = [F1; F2*1.25];
[t, ~, y_25_F2] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);
u = [F1; F2*1.5];
[t, ~, y_50_F2] = simulate_sde(t, xs, u, d, p, sigma_measurement, false, outputFolder);

du_10_F2 = (F2*1.1 - F2);
du_25_F2 = (F2*1.25 - F2);
du_50_F2 = (F2*1.5 - F2);

y_norm_10_F2 = (y_10_F2 - hs) / du_10_F2;
y_norm_25_F2 = (y_25_F2 - hs) / du_25_F2;
y_norm_50_F2 = (y_50_F2 - hs) / du_50_F2;

[fig_norm_F2, ax_norm_F2] = plot_step_responses(t, y_norm_10_F2, y_norm_25_F2, y_norm_50_F2, p, 'F_2');
sgtitle('Normalized Step Responses to F_2 Steps', 'FontSize', 14);
exportgraphics(fig_norm_F2, fullfile(outputFolder, 'normalized_steps_F2.pdf'), 'ContentType', 'vector');

%% Now trying system identification. with the noise it is better to do larger step. Showing results for 10% step (very bad) and 50% step

% 50% step
y_norm_F1 = (y_50_F1 - hs) / du_50_F1;
y_norm_F2 = (y_50_F2 - hs) / du_50_F2;

[fig_tf, ax_tf] = plot_transfer_identification(t, y_norm_F1, y_norm_F2);
outputFolder = fullfile('figures', 'problem_4', 'normalized_tf','sde', 'step_50'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(fig_tf, fullfile(outputFolder, 'normalized_tf.pdf'), 'ContentType', 'vector');

% 25% step
y_norm_F1 = (y_25_F1 - hs) / du_25_F1;
y_norm_F2 = (y_25_F2 - hs) / du_25_F2;

[fig_tf, ax_tf] = plot_transfer_identification(t, y_norm_F1, y_norm_F2);
outputFolder = fullfile('figures', 'problem_4', 'normalized_tf','sde', 'step_'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(fig_tf, fullfile(outputFolder, 'normalized_tf.pdf'), 'ContentType', 'vector');

% 10% step
y_norm_F1 = (y_10_F1 - hs) / du_10_F1;
y_norm_F2 = (y_10_F2 - hs) / du_10_F2;

[fig_tf, ax_tf] = plot_transfer_identification(t, y_norm_F1, y_norm_F2);
outputFolder = fullfile('figures', 'problem_4', 'normalized_tf','sde', 'step_10'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(fig_tf, fullfile(outputFolder, 'normalized_tf.pdf'), 'ContentType', 'vector');


[fig_tf, ax_tf] = plot_transfer_identification_combined(t, ...
    (y_10_F1 - hs)/du_10_F1, (y_25_F1 - hs)/du_25_F1, (y_50_F1 - hs)/du_50_F1, ...
    (y_10_F2 - hs)/du_10_F2, (y_25_F2 - hs)/du_25_F2, (y_50_F2 - hs)/du_50_F2, p);

outputFolder = fullfile('figures', 'problem_4', 'normalized_tf', 'sde', 'combined');
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(fig_tf, fullfile(outputFolder, 'normalized_tf_combined.pdf'), 'ContentType', 'vector');


%plot all combined. 4 plots in one figre. each has 3 step responses

% ------------------------------------------------------------
% Supporting residual function
% ------------------------------------------------------------
function r = FourTankSystemSteadyResidual(x, u, d, p)
r = FourTankSystemModified(0, x, u, d, p); % residual f(x)=0
end
