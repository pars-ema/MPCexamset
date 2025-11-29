%% ========================================================================
% PROBLEM 6 – Kalman Filters (Tasks 2, 3 & extra Task 4)
% Using ONLY the physical (Problem 5) linear model
%
% Task 2:
%   - Dynamic & static Kalman filters (non-augmented)
%   - Disturbances stochastic but NO step changes (Scenario A)
%   - Evaluation on linear plant with disturbances
%
% Task 3:
%   - Dynamic & static Kalman filters (AUGMENTED with [F3 F4] states)
%   - Disturbances stochastic AND contain step changes in F3 (Scenario B)
%   - Disturbance model for augmentation:
%       Random Walk  : d_{k+1} = d_k + w_d,k
%   - Evaluation on linear plant with disturbances
%
% Task 4 (extra):
%   - Step changes in inputs AND disturbances
%   - Nonlinear “truth” = SDE model (simulate_sde_2)
%   - Compare nonlinear vs linear vs KF (non-aug & aug)
%
% States are masses [g], outputs are heights [cm].
% ========================================================================

clear; clc; close all;

%% ========================================================================
% GLOBAL SETTINGS & OUTPUT FOLDER
% ========================================================================

Ts = 30;                % Sampling time [s]

outdir = fullfile('figures','problem_6');
if ~exist(outdir,'dir')
    mkdir(outdir);
end

%% ========================================================================
% 0) PHYSICAL PARAMETERS (FOR SDE + LINEARIZATION)
% ========================================================================

a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;        % [cm^2]
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327;% [cm^2]
g  = 981;                                                   % [cm/s^2]
rho = 1;                                                    % [g/cm^3]
gamma1 = 0.58;
gamma2 = 0.72;

p = [a1;a2;a3;a4;A1;A2;A3;A4;g;gamma1;gamma2;rho];

% Nominal inputs & disturbances (operating point)
u_s = [300;300];          % [cm^3/s] (F1, F2)
d_s = [250;250];          % [cm^3/s] (F3, F4)

% Steady-state masses xs for nonlinear SDE model
xs_guess = 500*ones(4,1);
opts = optimoptions('fsolve','Display','none');
xs = fsolve(@(x)FourTankSystemModified(0,x,u_s,d_s,p), xs_guess, opts);

% Steady-state heights (all 4 tanks)
h_s_abs = FourTankSystemSensor(xs,p);   % [4x1]

% We measure heights in tank 1 and 2
output_index = [1 2];
y_s_meas = h_s_abs(output_index);       % [2x1]

% Initial conditions (absolute masses)
x0 = xs;            % start at steady state
x0 = 0.8*xs;      % uncomment to start away from equilibrium

x0_dev = x0 - xs;   % deviation initial condition for linear model / KF

%% ========================================================================
% 1) CONTINUOUS LINEARIZATION & DISCRETIZATION (PHYSICAL MODEL)
% ========================================================================

nx = 4; nu = 2; nd = 2;

% Unpack parameters
a1=p(1); a2=p(2); a3=p(3); a4=p(4);
A1=p(5); A2=p(6); A3=p(7); A4=p(8);
g =p(9); gamma1=p(10); gamma2=p(11); rho=p(12);

% Torricelli constants (mass-based)
c1 = a1*sqrt(2*g/(rho*A1));
c2 = a2*sqrt(2*g/(rho*A2));
c3 = a3*sqrt(2*g/(rho*A3));
c4 = a4*sqrt(2*g/(rho*A4));

beta1 = c1/(2*sqrt(xs(1)));
beta2 = c2/(2*sqrt(xs(2)));
beta3 = c3/(2*sqrt(xs(3)));
beta4 = c4/(2*sqrt(xs(4)));

% Continuous-time A matrix (mass states)
A_c = rho * [
   -beta1,  0,       beta3,   0;
         0, -beta2,      0,  beta4;
         0,      0, -beta3,   0;
         0,      0,      0, -beta4
];

% Continuous-time B (pump flows F1,F2)
B_c = rho * [
    gamma1,        0;
         0,   gamma2;
         0, 1-gamma2;
   1-gamma1,       0
];

% Continuous-time disturbance input matrix E_c (F3,F4)
% Disturbances only affect tanks 3 and 4
E_c = [
    0, 0;
    0, 0;
    1, 0;
    0, 1
];

% Output matrix: masses -> h1,h2
C_c = [
    1/(rho*A1), 0, 0, 0;
    0, 1/(rho*A2), 0, 0
];

% Discretize [A_c, B_c, E_c] with ZOH using matrix exponential
M = [A_c, [B_c E_c];
     zeros(nu+nd, nx+nu+nd)];
Md = expm(M*Ts);

Ad = Md(1:nx, 1:nx);
BEd = Md(1:nx, nx+1:end);
Bd = BEd(:, 1:nu);
Ed = BEd(:, nu+1:end);     % [nx x nd]

Cd = C_c;                  % no direct discretization needed
Dd = zeros(size(Cd,1),nu);

% Use discrete model for KF design & linear simulation
A = Ad;
B = Bd;
C = Cd;
E = Ed;          % disturbance / process noise input on tanks 3 and 4

% Define G as the disturbance/noise input matrix (only tanks 3 and 4)
G = E;           % [4x2]

%% ========================================================================
% 2) TIME GRID + NOISE COVARIANCES
% ========================================================================

t  = 0:Ts:7200;
N  = numel(t);

% Inputs (Task 2 & 3 – constant at operating point)
u = repmat(u_s,1,N);    % [2xN]
u_dev = u - u_s;        % deviation input

% Measurement noise covariance R (heights h1,h2)
R = diag([4 4]);        % std = 2 cm on each level

% Process noise on disturbance channels w_k (2D: affects tanks 3 & 4)
Qw = diag([25 25]);     % tuning: variance on F3 and F4 "noise"

% Equivalent process noise on states: Qx = G Qw G'
Qx = G * Qw * G';

% Steady-state measured heights matrix
Y_s_meas_mat = repmat(y_s_meas,1,N);   % [2xN]

%% ========================================================================
% 3) TASK 2 – Non-augmented KFs on LINEAR PLANT with disturbances (Scenario A)
%     Disturbances stochastic, NO steps. Evaluation on linear model.
% ========================================================================

rng(0);   % reproducibility

% 3.1 Linear plant simulation with stochastic disturbances (no steps)
x_true_A = zeros(nx,N);
x_true_A(:,1) = x0_dev;     % deviation states

w_store = zeros(nd, N);
for k = 1:N-1
    w_k = mvnrnd(zeros(nd,1), Qw).';      % w_k ~ N(0,Qw)
    w_store(:,k) = w_k; 
    x_true_A(:,k+1) = A*x_true_A(:,k) ...
                    + B*u_dev(:,k) ...
                    + G*w_k;             % disturbances only on tanks 3 & 4
end

x_true_A_abs  = xs + x_true_A;                % absolute masses [g]
y_true_A_dev  = C*x_true_A;                  % deviation heights [cm]
y_true_A_abs  = Y_s_meas_mat + y_true_A_dev; % absolute heights [cm]

% Measurement noise (consistent with R)
y_meas_A_abs = y_true_A_abs ...
             + [sqrt(R(1,1))*randn(1,N);
                sqrt(R(2,2))*randn(1,N)];

% Deviation measurements for KF
y_meas_A_dev = y_meas_A_abs - Y_s_meas_mat;

% 3.2 Dynamic & static non-augmented KFs (designed for linear model)
KF_dyn_A = KF_dynamic(A,B,C,Qx,R,u_dev,y_meas_A_dev,x0_dev);
KF_sta_A = KF_static (A,B,C,Qx,R,u_dev,y_meas_A_dev,x0_dev);

% Deviation outputs from estimates
yhat_dyn_A_dev = C*KF_dyn_A.x_hat;
yhat_sta_A_dev = C*KF_sta_A.x_hat;

% Absolute outputs (heights)
yhat_dyn_A_abs = Y_s_meas_mat + yhat_dyn_A_dev;
yhat_sta_A_abs = Y_s_meas_mat + yhat_sta_A_dev;

% Absolute state estimates (masses)
x_dyn_A_abs = xs + KF_dyn_A.x_hat;
x_sta_A_abs = xs + KF_sta_A.x_hat;

% === Task 2 — Dynamic KF only ===
fig2_states_dyn = figure('Name','Task2_States_Dynamic');
sgtitle('States – Linear plant vs Dynamic KF');
colors = lines(2*nx);   % first nx colors = linear, next nx = KF
hold on; grid on;

for i = 1:nx
    % True linear plant (solid)
    plot(t, x_true_A_abs(i,:), '-',  'Color', colors(i,:), ...
         'LineWidth', 1.4, 'DisplayName', sprintf('Linear m_%d', i));

    % Dynamic KF estimate (dashed)
    plot(t, x_dyn_A_abs(i,:), '--', 'Color', colors(i+nx,:), ...
         'LineWidth', 1.7, 'DisplayName', sprintf('Dyn KF m_%d', i));
end

xlabel('Time [s]');
ylabel('Mass [g]');
xlim([0 2000]);
legend('Location','best');

exportgraphics(fig2_states_dyn, ...
    fullfile(outdir,'task2_states_masses_linear_dynamic.pdf'), ...
    'ContentType','vector');

%=== Task 2 — Static KF only ===
fig2_states_sta = figure('Name','Task2_States_Static');
sgtitle('States – Linear plant vs Static KF');

colors = lines(2*nx);   % first nx colors = linear, next nx = KF

hold on; grid on;

for i = 1:nx
    % True linear plant (solid)
    plot(t, x_true_A_abs(i,:), '-',  'Color', colors(i,:), ...
         'LineWidth', 1.4, 'DisplayName', sprintf('Linear m_%d', i));

    % Static KF estimate (dashed)
    plot(t, x_sta_A_abs(i,:), '--', 'Color', colors(i+nx,:), ...
         'LineWidth', 1.7, 'DisplayName', sprintf('Sta KF m_%d', i));
end

xlabel('Time [s]');
ylabel('Mass [g]');
xlim([0 2000]);
legend('Location','best');

exportgraphics(fig2_states_sta, ...
    fullfile(outdir,'task2_states_masses_linear_static.pdf'), ...
    'ContentType','vector');

% 3.4 TASK 2 – Outputs: linear plant vs measurements vs KFs
fig2_outputs = figure('Name','Task2_Outputs_Linear');
sgtitle('Outputs (heights) – Linear plant with stochastic disturbances');

% h1
subplot(2,1,1); hold on; grid on;
plot(t, y_true_A_abs(1,:),   'k','LineWidth',1.4,'DisplayName','True h_1 (linear)');
plot(t, y_meas_A_abs(1,:),   'k.','DisplayName','Measured h_1');
plot(t, yhat_dyn_A_abs(1,:), 'r','LineWidth',1.2,'DisplayName','Dyn KF h_1');
plot(t, yhat_sta_A_abs(1,:), 'b--','LineWidth',1.2,'DisplayName','Sta KF h_1');
ylabel('h_1 [cm]');
legend('Location','best');
xlim([0, 2000])

% h2
subplot(2,1,2); hold on; grid on;
plot(t, y_true_A_abs(2,:),   'k','LineWidth',1.4,'DisplayName','True h_2 (linear)');
plot(t, y_meas_A_abs(2,:),   'k.','DisplayName','Measured h_2');
plot(t, yhat_dyn_A_abs(2,:), 'r','LineWidth',1.2,'DisplayName','Dyn KF h_2');
plot(t, yhat_sta_A_abs(2,:), 'b--','LineWidth',1.2,'DisplayName','Sta KF h_2');
ylabel('h_2 [cm]'); xlabel('Time [s]');
legend('Location','best');
xlim([0, 2000])

set(fig2_outputs,'Position',[100 100 800 500]);

exportgraphics(fig2_outputs, ...
    fullfile(outdir,'task2_outputs_heights_linear.pdf'), ...
    'ContentType','vector');

%% ========================================================================
% 4) TASK 3 – Non-aug & augmented KFs on LINEAR PLANT (Scenario B)
%     Disturbances stochastic WITH step in F3
% ========================================================================

rng(1);   % different seed

% 4.1 Linear plant simulation with disturbance step in F3 (channel 1 of w)
x_true_B = zeros(nx,N);
x_true_B(:,1) = x0_dev;

mu_w = zeros(nd,N);         % mean of w_k
mu_step = 15;               % step magnitude in disturbance channel 1 [cm^3/s]
k_step = floor(3000/Ts);
mu_w(1, k_step:end) = mu_step;

for k = 1:N-1
    w_k = mu_w(:,k) + mvnrnd(zeros(nd,1), Qw).';  % w_k ~ N(mu_w, Qw)
    x_true_B(:,k+1) = A*x_true_B(:,k) ...
                    + B*u_dev(:,k) ...
                    + G*w_k;
end

x_true_B_abs  = xs + x_true_B;
y_true_B_dev  = C*x_true_B;
y_true_B_abs  = Y_s_meas_mat + y_true_B_dev;

% Measurement noise
y_meas_B_abs = y_true_B_abs ...
             + [sqrt(R(1,1))*randn(1,N);
                sqrt(R(2,2))*randn(1,N)];
y_meas_B_dev = y_meas_B_abs - Y_s_meas_mat;

% 4.2 Non-augmented KFs (do NOT model disturbance state)
KF_dyn_B_nonaug = KF_dynamic(A,B,C,Qx,R,u_dev,y_meas_B_dev,x0_dev);
KF_sta_B_nonaug = KF_static (A,B,C,Qx,R,u_dev,y_meas_B_dev,x0_dev);

yhat_dyn_B_nonaug_dev = C*KF_dyn_B_nonaug.x_hat;
yhat_sta_B_nonaug_dev = C*KF_sta_B_nonaug.x_hat;

yhat_dyn_B_nonaug_abs = Y_s_meas_mat + yhat_dyn_B_nonaug_dev;
yhat_sta_B_nonaug_abs = Y_s_meas_mat + yhat_sta_B_nonaug_dev;

x_dyn_B_nonaug_abs = xs + KF_dyn_B_nonaug.x_hat;
x_sta_B_nonaug_abs = xs + KF_sta_B_nonaug.x_hat;

% 4.3 Augmented dynamic & static KFs (estimate disturbance states)
Qd_RW = 1 * eye(nd);      % disturbance covariance for random-walk model

KF_dyn_B_aug = KF_aug_dynamic(A,B,G,C,Qx,Qd_RW,R,u_dev,y_meas_B_dev,x0_dev);
KF_sta_B_aug = KF_aug_static (A,B,G,C,Qx,Qd_RW,R,u_dev,y_meas_B_dev,x0_dev);

x_dyn_B_aug_dev = KF_dyn_B_aug.x_hat(1:nx,:);        % deviation states
d_dyn_B_aug     = KF_dyn_B_aug.x_hat(nx+1:end,:);    % [ΔF3; ΔF4]

x_sta_B_aug_dev = KF_sta_B_aug.x_hat(1:nx,:);
d_sta_B_aug     = KF_sta_B_aug.x_hat(nx+1:end,:);

x_dyn_B_aug_abs = xs + x_dyn_B_aug_dev;
x_sta_B_aug_abs = xs + x_sta_B_aug_dev;

yhat_dyn_B_aug_dev = C*x_dyn_B_aug_dev;
yhat_sta_B_aug_dev = C*x_sta_B_aug_dev;

yhat_dyn_B_aug_abs = Y_s_meas_mat + yhat_dyn_B_aug_dev;
yhat_sta_B_aug_abs = Y_s_meas_mat + yhat_sta_B_aug_dev;

% 4.4 TASK 3 – Outputs (linear plant): non-aug vs augmented
fig3_outputs = figure('Name','Task3_Outputs_Linear');
sgtitle('Outputs – Linear plant with disturbance step in F_3');

% h1
subplot(2,1,1); hold on; grid on;
plot(t, y_true_B_abs(1,:),          'k','LineWidth',1.4,'DisplayName','True h_1 (linear)');
plot(t, y_meas_B_abs(1,:),          'k.','DisplayName','Measured h_1');
plot(t, yhat_dyn_B_nonaug_abs(1,:),'r','LineWidth',1.2,'DisplayName','Dyn non-aug');
plot(t, yhat_dyn_B_aug_abs(1,:),   'm','LineWidth',1.4,'DisplayName','Dyn aug');
plot(t, yhat_sta_B_aug_abs(1,:),   'g--','LineWidth',1.4,'DisplayName','Sta aug');
legend('Location','best');
xline(t(k_step), 'k--', 'Step in F_3 disturbance', 'LabelVerticalAlignment','bottom', FontSize=10, Label='Step in F3');
ylabel('h_1 [cm]');


% h2
subplot(2,1,2); hold on; grid on;
plot(t, y_true_B_abs(2,:),          'k','LineWidth',1.4,'DisplayName','True h_2 (linear)');
plot(t, y_meas_B_abs(2,:),          'k.','DisplayName','Measured h_2');
plot(t, yhat_dyn_B_nonaug_abs(2,:),'r','LineWidth',1.2,'DisplayName','Dyn non-aug');
plot(t, yhat_sta_B_nonaug_abs(2,:),'r:','LineWidth',1.2,'DisplayName','Sta non-aug');
plot(t, yhat_dyn_B_aug_abs(2,:),   'm','LineWidth',1.4,'DisplayName','Dyn aug');
plot(t, yhat_sta_B_aug_abs(2,:),   'g--','LineWidth',1.4,'DisplayName','Sta aug');
xline(t(k_step), 'k--', 'Step in F_3 disturbance', 'LabelVerticalAlignment','bottom', FontSize=10, Label='Step in F3');
ylabel('h_2 [cm]'); xlabel('Time [s]');
legend('Location','best');

set(fig3_outputs,'Position',[100 100 800 500]);

exportgraphics(fig3_outputs, ...
    fullfile(outdir,'task3_outputs_heights_linear.pdf'), ...
    'ContentType','vector');

%% 4.5 TASK 3 – States in tank layout [3 4; 1 2]

fig3_states = figure('Name','Task3_States_Linear');
sgtitle('States – Linear plant vs KFs (Tank Layout)');

% Tank index ordering corresponding to physical layout
tank_order = [3 4 1 2];

for o = 1:4
    i = tank_order(o);   % actual state index (mass index)
    subplot(2,2,o); hold on; grid on;

    % True linear plant
    plot(t, x_true_B_abs(i,:), 'k', 'LineWidth',1.4, ...
        'DisplayName', sprintf('True m_%d (linear)',i));

    % Non-augmented KF
    plot(t, x_dyn_B_nonaug_abs(i,:), 'r',  'LineWidth',1.2, ...
        'DisplayName','Dyn non-aug');

    % Augmented KF
    plot(t, x_dyn_B_aug_abs(i,:), 'm',   'LineWidth',1.3, ...
        'DisplayName','Dyn aug');
    plot(t, x_sta_B_aug_abs(i,:), 'g--', 'LineWidth',1.3, ...
        'DisplayName','Sta aug');

    ylabel(sprintf('m_%d [g]', i));
    title(sprintf('Tank %d', i));

    % Legend for each subplot
    legend('Location','best');
    xlim([0 5000])
end

xlabel('Time [s]');
set(fig3_states,'Position',[100 100 900 700]);

exportgraphics(fig3_states, ...
    fullfile(outdir,'task3_states_masses_linear_layout.pdf'), ...
    'ContentType','vector');

%% 4.6 TASK 3 – Disturbance estimates (true, mean, and KF estimates)

fig3_dist = figure('Name','Task3_Disturbances_Linear');
sgtitle('Disturbances – Mean of the disturbances and their estimates');

% Extract actual disturbances injected in the linear plant
d_true_F3 = w_store(1,:);     % actual injected disturbance for F3
d_true_F4 = w_store(2,:);     % actual injected disturbance for F4

% Extract mean disturbances (applied step)
mean_F3 = mu_w(1,:);          % mean bias for F3
mean_F4 = mu_w(2,:);          % mean bias for F4

% --- F3 subplot ---
subplot(2,1,1); hold on; grid on;
plot(t, mean_F3,            'k',  'LineWidth',1.8, 'DisplayName','Mean F_3');
plot(t, d_dyn_B_aug(1,:),   'm',  'LineWidth',1.4, 'DisplayName','Dyn KF (augmented)');
plot(t, d_sta_B_aug(1,:),   'g--','LineWidth',1.4, 'DisplayName','Sta KF (augmented)');
ylabel('\Delta F_3 [cm^3/s]');
title('Disturbance F_3 (Tank 3)');
legend('Location','best');
xlim([0 t(end)]);
ylim([-20 20]);

% --- F4 subplot ---
subplot(2,1,2); hold on; grid on;
plot(t, mean_F4,            'k',  'LineWidth',1.8, 'DisplayName','Mean F_4');
plot(t, d_dyn_B_aug(2,:),   'm',  'LineWidth',1.4, 'DisplayName','Dyn KF (augmented)');
plot(t, d_sta_B_aug(2,:),   'g--','LineWidth',1.4, 'DisplayName','Sta KF (augmented)');
ylabel('\Delta F_4 [cm^3/s]');
xlabel('Time [s]');
title('Disturbance F_4 (Tank 4)');
legend('Location','best');
xlim([0 t(end)]);
ylim([-20 20]);

set(fig3_dist,'Position',[100 100 850 450]);

exportgraphics(fig3_dist, ...
    fullfile(outdir,'task3_disturbances_linear_mean.pdf'), ...
    'ContentType','vector');

%% ========================================================================
% 5) TASK 4 – Step Changes in Inputs AND Disturbances (extra)
%     Nonlinear SDE truth vs linear vs KF (non-aug & aug)
% ========================================================================

disp('Running Task 4: step changes in inputs and disturbances (SDE truth)...');

% 5.1 Step-changing inputs u1, u2
u_var = repmat(u_s,1,N);

x0 = [0;0;0;0];
x0_dev = x0 - xs;

k1 = floor(1000/Ts);
k2 = floor(2500/Ts);
k3 = floor(4500/Ts);
k4 = floor(5000/Ts);

% u1 steps
u_var(1, k1:end) = u_s(1) + 50;     % +50 at t≈1000s
u_var(1, k2:end) = u_s(1) - 30;     % -30 at t≈2500s

% u2 steps
u_var(2, k2:end) = u_s(2) + 40;     % +40 at t≈2500s
u_var(2, k3:end) = u_s(2) - 20;     % -20 at t≈4500s

u_var_dev = u_var - u_s;

% 5.2 Step-changing disturbances F3,F4 (absolute)
d_var = repmat(d_s,1,N);
d_var(1, k1:end) = d_s(1) + 30;  
d_var(1, k3:end) = d_s(1) - 20;  
d_var(2, k2:end) = d_s(2) + 50;
d_var(2, k4:end) = d_s(2) + 10;

% 5.3 Plot u and d steps
fig4_ud = figure('Name','Task4_InputDisturbanceSteps');
sgtitle('Step Changes in Inputs and Disturbances');

subplot(2,1,1); hold on; grid on;
plot(t,u_var(1,:),'b','LineWidth',1.4);
plot(t,u_var(2,:),'r','LineWidth',1.4);
ylabel('u_1, u_2 [cm^3/s]');
legend('u_1','u_2');
title('Inputs');

subplot(2,1,2); hold on; grid on;
plot(t,d_var(1,:),'b','LineWidth',1.4);
plot(t,d_var(2,:),'r','LineWidth',1.4);
ylabel('F_3, F_4 [cm^3/s]');
legend('F_3','F_4');
xlabel('Time [s]');
title('Disturbances');

set(fig4_ud,'Position',[100 100 800 500]);

exportgraphics(fig4_ud, ...
    fullfile(outdir,'task4_u_d_steps.pdf'), ...
    'ContentType','vector');

% 5.4 Nonlinear SDE simulation (truth)
sigma_meas = [2;2;2;2];           % [cm], for SDE simulator

[x_true_var, y_true_var_abs, y_meas_var_abs, y_meas_var_dev] = simulate_sde_2( ...
    t, x0, u_var, d_var, p, sigma_meas, output_index, xs);

y_true_var_12 = y_true_var_abs(output_index,:);

% 5.5 Linear model simulation (with disturbance input)
x_lin_var = zeros(nx,N);
x_lin_var(:,1) = x0_dev;

for k = 1:N-1
    x_lin_var(:,k+1) = A*x_lin_var(:,k) ...
                     + B*u_var_dev(:,k) ...
                     + E*(d_var(:,k)-d_s);  % explicit disturbance deviations
end

x_lin_var_abs = xs + x_lin_var;
y_lin_var_abs = Y_s_meas_mat + C*x_lin_var;

% 5.6 Kalman Filters for Task 4 (dynamic only)
KF_dyn_var_nonaug = KF_dynamic(A,B,C,Qx,R,u_var_dev,y_meas_var_dev,x0_dev);

KF_dyn_var_aug = KF_aug_dynamic(A,B,G,C,Qx,Qd_RW,R,u_var_dev, ...
                                y_meas_var_dev,x0_dev);

x_dyn_var_nonaug_abs = xs + KF_dyn_var_nonaug.x_hat;
x_dyn_var_aug_abs    = xs + KF_dyn_var_aug.x_hat(1:nx,:);

y_dyn_var_nonaug_abs = Y_s_meas_mat + C*KF_dyn_var_nonaug.x_hat;
y_dyn_var_aug_abs    = Y_s_meas_mat + C*KF_dyn_var_aug.x_hat(1:nx,:);

d_dyn_var = KF_dyn_var_aug.x_hat(nx+1:end,:);   % disturbance estimates

% 5.7 TASK 4 – Outputs
fig4_outputs = figure('Name','Task4_Outputs');
sgtitle('Outputs with Step Changes – SDE vs linear vs KFs');

subplot(2,1,1); hold on; grid on;
plot(t,y_true_var_12(1,:),'k','LineWidth',1.4,'DisplayName','True h_1 (SDE)');
plot(t,y_meas_var_abs(1,:),'k.','DisplayName','Measured h_1');
plot(t,y_lin_var_abs(1,:),'c--','LineWidth',1.1,'DisplayName','Linear h_1');
plot(t,y_dyn_var_nonaug_abs(1,:),'r','LineWidth',1.2,'DisplayName','KF non-aug');
plot(t,y_dyn_var_aug_abs(1,:),'m','LineWidth',1.4,'DisplayName','KF aug');
ylabel('h_1 [cm]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(t,y_true_var_12(2,:),'k','LineWidth',1.4,'DisplayName','True h_2 (SDE)');
plot(t,y_meas_var_abs(2,:),'k.','DisplayName','Measured h_2');
plot(t,y_lin_var_abs(2,:),'c--','LineWidth',1.1,'DisplayName','Linear h_2');
plot(t,y_dyn_var_nonaug_abs(2,:),'r','LineWidth',1.2,'DisplayName','KF non-aug');
plot(t,y_dyn_var_aug_abs(2,:),'m','LineWidth',1.4,'DisplayName','KF aug');
ylabel('h_2 [cm]'); xlabel('Time [s]');
legend('Location','best');

set(fig4_outputs,'Position',[100 100 800 500]);

exportgraphics(fig4_outputs, ...
    fullfile(outdir,'task4_outputs_heights.pdf'), ...
    'ContentType','vector');

% 5.8 TASK 4 – States
fig4_states = figure('Name','Task4_States');
sgtitle('State Estimation with Steps – SDE vs linear vs KFs');

for i = 1:nx
    subplot(nx,1,i); hold on; grid on;
    plot(t,x_true_var(i,:),'k','LineWidth',1.4,'DisplayName','True SDE');
    plot(t,x_lin_var_abs(i,:),'c--','LineWidth',1.1,'DisplayName','Linear model');
    plot(t,x_dyn_var_nonaug_abs(i,:),'r','LineWidth',1.2,'DisplayName','KF non-aug');
    plot(t,x_dyn_var_aug_abs(i,:),'m','LineWidth',1.3,'DisplayName','KF aug');
    ylabel(sprintf('m_%d [g]',i));
    if i==1
        legend('Location','best');
    end
end
xlabel('Time [s]');
set(fig4_states,'Position',[100 100 800 600]);

exportgraphics(fig4_states, ...
    fullfile(outdir,'task4_states_masses.pdf'), ...
    'ContentType','vector');

% 5.9 TASK 4 – Disturbance estimates
fig4_dist = figure('Name','Task4_Disturbances');
sgtitle('Estimated Disturbances');

subplot(2,1,1); hold on; grid on;
plot(t,d_var(1,:) - d_s(1),'k','LineWidth',1.4,'DisplayName','True \DeltaF_3');
plot(t,d_dyn_var(1,:),'m','LineWidth',1.4,'DisplayName','KF \DeltaF_3');
ylabel('\DeltaF_3 [cm^3/s]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(t,d_var(2,:) - d_s(2),'k','LineWidth',1.4,'DisplayName','True \DeltaF_4');
plot(t,d_dyn_var(2,:),'m','LineWidth',1.4,'DisplayName','KF \DeltaF_4');
ylabel('\DeltaF_4 [cm^3/s]'); xlabel('Time [s]');
legend('Location','best');

set(fig4_dist,'Position',[100 100 800 500]);

exportgraphics(fig4_dist, ...
    fullfile(outdir,'task4_disturbances.pdf'), ...
    'ContentType','vector');

disp('Problem 6 simulations finished. PDFs saved in figures/problem_6.');

%% ========================================================================
% ===================== KALMAN FILTER FUNCTIONS ===========================
% ========================================================================

function KF = KF_dynamic(A,B,C,Q,R,u,y,x0)
    % Dynamic (time-varying) Kalman filter for deviation model:
    %   x_{k+1} = A x_k + B u_k + w_k,   w ~ N(0,Q)
    %   y_k     = C x_k + v_k,           v ~ N(0,R)
    nx = size(A,1);
    N  = size(y,2);

    x = zeros(nx,N);
    x(:,1) = x0;
    P = 1e-2*eye(nx);

    for k = 1:N-1
        % Measurement update
        S = C*P*C' + R;
        K = P*C'/S;
        innov = y(:,k) - C*x(:,k);
        x(:,k) = x(:,k) + K*innov;
        P = (eye(nx)-K*C)*P;

        % Time update
        x(:,k+1) = A*x(:,k) + B*u(:,k);
        P = A*P*A' + Q;
    end

    KF.x_hat = x;
end

function KF = KF_static(A,B,C,Q,R,u,y,x0)
    % STATIC (steady-state) Kalman filter:
    % Gain K is constant, obtained from DARE.
    nx = size(A,1);
    N  = size(y,2);

    [P_ss,~,~] = dare(A',C',Q,R);
    P_ss = 0.5*(P_ss+P_ss');    % enforce symmetry

    Re = C*P_ss*C' + R;
    K  = P_ss*C'/Re;

    x = zeros(nx,N);
    x(:,1) = x0;

    for k = 1:N-1
        innov = y(:,k) - C*x(:,k);
        x(:,k) = x(:,k) + K*innov;
        x(:,k+1) = A*x(:,k) + B*u(:,k);
    end

    KF.x_hat = x;
    KF.K     = K;
    KF.P     = P_ss;
end

function KF = KF_aug_dynamic(A,B,G,C,Qx,Qd,R,u,y,x0_dev)
    % Builds augmented model (x + disturbance states) and runs dynamic KF.
    % Disturbance model: Random Walk
    %
    % Original model:
    %   x_{k+1} = A x_k + B u_k + G d_k + w_x,k
    % Augmented model:
    %   [x_{k+1}] = [A  G][x_k] + [B]u_k + [w_x,k]
    %   [d_{k+1}]   [0  I][d_k]   [0]      [w_d,k]
    % y_k = [C  0][x_k;d_k] + v_k

    nx = size(A,1);
    nd = size(G,2);
    ny = size(C,1);
    nu = size(B,2);

    Aa = [A, G;
          zeros(nd,nx), eye(nd)];
    Ba = [B;
          zeros(nd,nu)];
    Ca = [C, zeros(ny,nd)];
    Qa = blkdiag(Qx,Qd);

    x0_aug = [x0_dev; zeros(nd,1)];
    KF = KF_dynamic(Aa,Ba,Ca,Qa,R,u,y,x0_aug);
end

function KF = KF_aug_static(A,B,G,C,Qx,Qd,R,u,y,x0_dev)
    % Builds augmented model (x + disturbance states) and runs static KF.
    % Disturbance model: Random Walk
    nx = size(A,1);
    nd = size(G,2);
    ny = size(C,1);
    nu = size(B,2);

    Aa = [A, G;
          zeros(nd,nx), eye(nd)];
    Ba = [B;
          zeros(nd,nu)];
    Ca = [C, zeros(ny,nd)];
    Qa = blkdiag(Qx,Qd);

    x0_aug = [x0_dev; zeros(nd,1)];
    KF = KF_static(Aa,Ba,Ca,Qa,R,u,y,x0_aug);
end


