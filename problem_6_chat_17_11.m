%% ===============================================================
% Problem 6 – Kalman Filters for Linearized Four-Tank System
% Using OUTPUT-BIAS disturbance augmentation
%
% States are MASS deviations (x = m - m_s).
%
% Filters:
%   • Dynamic KF (non-augmented)
%   • Static  KF (non-augmented)
%   • Dynamic augmented KF (x + bias b)
%   • Static  augmented KF (x + bias b)
%
% Scenarios:
%   • A: stochastic disturbances, NO step
%   • B: stochastic disturbances, WITH 10% step in F3
%
% Plant model  = nonlinear SDE (simulate_sde_2)
% Observer     = discrete linear model (Ad,Bd,Cd) from linearization
%
% Requires:
%   FourTankSystemModified.m
%   simulate_sde_2.m
%   discrete_model_from_linearization.mat  (Ad, Bd, Cd, Dd)
% ===============================================================
clc; clear; close all;

%% ===============================================================
% 1) Parameters and steady state of nonlinear model
% ===============================================================
a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327;
g  = 981; 
rho = 1; 
gamma1 = 0.58; 
gamma2 = 0.72;

p = [a1;a2;a3;a4;A1;A2;A3;A4;g;gamma1;gamma2;rho];

u_s = [300;300];    % steady inputs [F1;F2]
d_s = [250;250];    % steady disturbances [F3;F4]

% Steady-state masses xs from nonlinear model
xs0  = 500*ones(4,1);
opts = optimoptions('fsolve','Display','none');
xs   = fsolve(@(x)FourTankSystemModified(0,x,u_s,d_s,p),xs0,opts);

% Corresponding heights
hs = [ xs(1)/(rho*A1);
       xs(2)/(rho*A2);
       xs(3)/(rho*A3);
       xs(4)/(rho*A4) ];
ys = hs(1:2);   % measured outputs (heights of tank 1 & 2)

disp('Steady-state masses xs:');
disp(xs.');
disp('Steady-state heights hs:');
disp(hs.');

%% ===============================================================
% 2) Load discrete linear model from linearization
% ===============================================================
load('discrete_model_from_linearization.mat','Ad','Bd','Cd','Dd');

[nx,nu] = size(Bd);   % expected 4 x 2
ny      = size(Cd,1); % expected 2

disp('Loaded discrete linear model:');
Ad, Bd, Cd

Ts = 30;          % sampling time [s]
t  = 0:Ts:7200;   % simulation horizon
N  = numel(t);

% Inputs (absolute, constant at u_s; deviation du defined later)
u = repmat(u_s,1,N);

%% ===============================================================
% 3) Noise covariances for KFs
% ===============================================================
% Here Qx is directly the covariance on the state deviation x.
% (Equivalent to using G = I and Q_w = Qx in the theory.)
Qx = 1e-4 * eye(nx);        % process noise covariance on x
R  = diag([2;2].^2);        % measurement noise covariance on y (h1,h2)

% Output bias random walk covariance (for augmented filters)
Qb_small = 1e-3 * eye(ny);  % scenario A: small drift (if used)
Qb_big   = 1e-2  * eye(ny);  % scenario B: allow fast bias change

%% Helper: mass->height conversion
mass2height = @(x) [ x(1,:)/(rho*A1);
                     x(2,:)/(rho*A2);
                     x(3,:)/(rho*A3);
                     x(4,:)/(rho*A4) ];

%% ===============================================================
% 4) SCENARIO A – stochastic disturbances, NO step in F3
% ===============================================================
disp('--- SCENARIO A: stochastic disturbances, NO step in F3 ---');

% Constant disturbances at steady-state
dA = repmat(d_s,1,N);

% ---- Nonlinear SDE simulation (true plant) ----
% Assumed signature: [x_sde, y_true, y_meas] = simulate_sde_2(...)
[x_sde_A, y_true_A_full, y_meas_A_full] = simulate_sde_2( ...
    t, xs, u, dA, p, [2;2;2;2], false, []);

% True states and heights (absolute)
x_true_A = x_sde_A;              % 4 x N masses (absolute)
h_true_A = y_true_A_full;        % assume 4 x N heights from SDE

% Keep only measured outputs (h1,h2)
y_true_A = y_true_A_full(1:2,:); % 2 x N
y_meas_A = y_meas_A_full(1:2,:); % 2 x N

% Deviation variables for KF
dy_meas_A = y_meas_A - ys;       % measurement deviation 2 x N
du        = u - u(:,1);          % input deviation 2 x N

% ---- Linear discrete model simulation of deviations ----
[x_lin_A_dev, y_lin_A] = simulate_linear_ss(Ad,Bd,Cd,du,ys,N);
x_lin_A = x_lin_A_dev + xs;      % absolute masses from linear model
h_lin_A = mass2height(x_lin_A);  % absolute heights from linear model

% ---- Dynamic & Static Kalman filters (non-augmented) ----
resA = runKF_nonaug(Ad,Bd,Cd,Qx,R,dy_meas_A,du,ys);

%% ===============================================================
% 5) Plots – Scenario A
% ===============================================================

% 5.1 Outputs h1,h2: true vs measured vs linear vs KFs
figure;
subplot(2,1,1); hold on; grid on;
plot(t, y_true_A(1,:),   'k',  'LineWidth', 1.5);
plot(t, y_meas_A(1,:),   'bo', 'MarkerSize', 3);
plot(t, y_lin_A(1,:),    'c--','LineWidth', 1.3);
plot(t, resA.y_dyn(1,:), 'r',  'LineWidth', 1.3);
plot(t, resA.y_sta(1,:), 'g--','LineWidth', 1.3);
xlabel('Time [s]'); ylabel('h_1 [cm]');
title('Scenario A – Output h_1');
legend('True (SDE)','Measured','Linear model','Dyn KF','Static KF','Location','Best');

subplot(2,1,2); hold on; grid on;
plot(t, y_true_A(2,:),   'k',  'LineWidth', 1.5);
plot(t, y_meas_A(2,:),   'bo', 'MarkerSize', 3);
plot(t, y_lin_A(2,:),    'c--','LineWidth', 1.3);
plot(t, resA.y_dyn(2,:), 'r',  'LineWidth', 1.3);
plot(t, resA.y_sta(2,:), 'g--','LineWidth', 1.3);
xlabel('Time [s]'); ylabel('h_2 [cm]');
title('Scenario A – Output h_2');

% 5.2 State estimates vs true states (heights) – dynamic & static KF
h_dyn_A = mass2height(resA.x_dyn + xs);   % KF absolute heights
h_sta_A = mass2height(resA.x_sta + xs);   % KF absolute heights

figure;
for i = 1:4
    subplot(2,2,i); hold on; grid on;
    plot(t, h_true_A(i,:), 'k',  'LineWidth', 1.5); % true heights
    plot(t, h_lin_A(i,:),  'c--','LineWidth', 1.2); % linear model heights
    plot(t, h_dyn_A(i,:),  'r',  'LineWidth', 1.2); % dyn KF heights
    plot(t, h_sta_A(i,:),  'g--','LineWidth', 1.2); % static KF heights
    xlabel('Time [s]');
    ylabel(sprintf('h_%d [cm]',i));
    title(sprintf('Scenario A – State h_%d',i));
    if i == 1
        legend('True (SDE)','Linear model','Dyn KF','Static KF','Location','Best');
    end
end

% 5.3 State estimates vs true states (masses) – dynamic & static KF
m_true_A = x_true_A;                 % true masses (absolute)
m_dyn_A  = resA.x_dyn + xs;          % KF masses (absolute)
m_sta_A  = resA.x_sta + xs;          % KF masses (absolute)
m_lin_A  = x_lin_A;                  % linear model masses (absolute)

figure;
for i = 1:4
    subplot(2,2,i); hold on; grid on;
    plot(t, m_true_A(i,:), 'k',  'LineWidth', 1.5); % true masses
    plot(t, m_lin_A(i,:),  'c--','LineWidth', 1.2); % linear model
    plot(t, m_dyn_A(i,:),  'r',  'LineWidth', 1.2); % dyn KF
    plot(t, m_sta_A(i,:),  'g--','LineWidth', 1.2); % static KF
    xlabel('Time [s]');
    ylabel(sprintf('m_%d [g]',i));
    title(sprintf('Scenario A – State m_%d',i));
    if i == 1
        legend('True (SDE)','Linear model','Dyn KF','Static KF','Location','Best');
    end
end

%% ===============================================================
% 6) SCENARIO B – stochastic disturbances + 10% step in F3
% ===============================================================
disp('--- SCENARIO B: stochastic disturbances, 10% step in F3 ---');

dB = repmat(d_s,1,N);
k_step = find(t >= 3000, 1);
dB(1,k_step:end) = 1.10 * d_s(1);   % 10% step in F3 at t = 3000 s

% Nonlinear SDE simulation (true plant)
[x_sde_B, y_true_B_full, y_meas_B_full] = simulate_sde_2( ...
    t, xs, u, dB, p, [4;4;4;4], false, []);

x_true_B = x_sde_B;                      % true masses (absolute)
h_true_B = mass2height(x_true_B);        % true heights (computed)

y_true_B = y_true_B_full(1:2,:);
y_meas_B = y_meas_B_full(1:2,:);

dy_meas_B = y_meas_B - ys;   % deviation
% du is same as before (inputs constant around u_s)

% Linear discrete model response (deviation)
[x_lin_B_dev, y_lin_B] = simulate_linear_ss(Ad,Bd,Cd,du,ys,N);
x_lin_B = x_lin_B_dev + xs;          % absolute masses
h_lin_B = mass2height(x_lin_B);      % absolute heights

% Non-augmented KFs
resB_nonaug = runKF_nonaug(Ad,Bd,Cd,Qx,R,dy_meas_B,du,ys);

% Augmented KFs (x + bias b)
resB_aug    = runKF_aug(Ad,Bd,Cd,Qx,Qb_big,R,dy_meas_B,du,ys);

%% ===============================================================
% 7) Plots – Scenario B
% ===============================================================

% 7.1 Outputs h1,h2: compare everything
fig = figure('Units','normalized','Position',[0 0 0.75 0.75]);
% (your plotting code...)

subplot(2,1,1); hold on; grid on;
plot(t, y_true_B(1,:),          'k',  'LineWidth', 1.5);
plot(t, y_meas_B(1,:),          'bo', 'MarkerSize', 3);
plot(t, y_lin_B(1,:),           'c--','LineWidth', 1.2);
plot(t, resB_nonaug.y_dyn(1,:)+0.1, 'r:', 'LineWidth', 1.2);
plot(t, resB_nonaug.y_sta(1,:), 'g:', 'LineWidth', 1.2);
plot(t, resB_aug.y_dyn(1,:),    'm',  'LineWidth', 1.4);
plot(t, resB_aug.y_sta(1,:),    'b--','LineWidth', 1.4);
xline(t(k_step),'k--','Step in F_3');
xlabel('Time [s]'); ylabel('h_1 [cm]');
title('Scenario B – Output h_1 (10% step in F_3)');
legend('True (SDE)','Measured','Linear model', ...
       'Dyn KF (non-aug)','Static KF (non-aug)', ...
       'Dyn KF (aug)','Static KF (aug)','Location','Best');

subplot(2,1,2); hold on; grid on;
plot(t, y_true_B(2,:),          'k',  'LineWidth', 1.5);
plot(t, y_meas_B(2,:),          'bo', 'MarkerSize', 3);
plot(t, y_lin_B(2,:),           'c--','LineWidth', 1.2);
plot(t, resB_nonaug.y_dyn(2,:), 'r:', 'LineWidth', 1.2);
plot(t, resB_nonaug.y_sta(2,:), 'g:', 'LineWidth', 1.2);
plot(t, resB_aug.y_dyn(2,:),    'm',  'LineWidth', 1.4);
plot(t, resB_aug.y_sta(2,:),    'b--','LineWidth', 1.4);
xline(t(k_step),'k--','Step in F_3');
xlabel('Time [s]'); ylabel('h_2 [cm]');
title('Scenario B – Output h_2 (10% step in F_3)');
exportgraphics(gcf, 'small_q.pdf', 'ContentType', 'vector');

% 7.2 Bias estimates b1,b2
figure; hold on; grid on;
plot(t, resB_aug.b_dyn(1,:), 'r',  'LineWidth', 1.5);
plot(t, resB_aug.b_sta(1,:), 'g--','LineWidth', 1.5);
xline(t(k_step),'k--','Step in F_3');
xlabel('Time [s]'); ylabel('b_1 [cm]');
title('Scenario B – Estimated output bias b_1');
legend('Dyn Aug KF','Static Aug KF','Step in F_3','Location','Best');

figure; hold on; grid on;
plot(t, resB_aug.b_dyn(2,:), 'r',  'LineWidth', 1.5);
plot(t, resB_aug.b_sta(2,:), 'g--','LineWidth', 1.5);
xline(t(k_step),'k--','Step in F_3');
xlabel('Time [s]'); ylabel('b_2 [cm]');
title('Scenario B – Estimated output bias b_2');
legend('Dyn Aug KF','Static Aug KF','Step in F_3','Location','Best');

% 7.3 State estimates (heights) – augmented KFs
h_dyn_B_aug = mass2height(resB_aug.x_dyn + xs);
h_sta_B_aug = mass2height(resB_aug.x_sta + xs);

figure;
for i = 1:4
    subplot(2,2,i); hold on; grid on;
    plot(t, h_true_B(i,:),      'k',  'LineWidth', 1.5); % true heights
    plot(t, h_lin_B(i,:),       'c--','LineWidth', 1.2); % linear model
    plot(t, h_dyn_B_aug(i,:),   'm',  'LineWidth', 1.2); % dyn aug KF
    plot(t, h_sta_B_aug(i,:),   'b--','LineWidth', 1.2); % static aug KF
    xline(t(k_step),'k--','Step in F_3');
    xlabel('Time [s]');
    ylabel(sprintf('h_%d [cm]',i));
    title(sprintf('Scenario B – State h_%d',i));
    if i == 1
        legend('True (SDE)','Linear model', ...
               'Dyn Aug KF','Static Aug KF','Step in F_3','Location','Best');
    end
end

% 7.4 State estimates (masses) – augmented KFs
m_true_B    = x_true_B;              % true masses (absolute)
m_dyn_B_aug = resB_aug.x_dyn + xs;   % dyn aug KF masses (absolute)
m_sta_B_aug = resB_aug.x_sta + xs;   % static aug KF masses (absolute)
m_lin_B     = x_lin_B;               % linear model masses (absolute)

figure;
for i = 1:4
    subplot(2,2,i); hold on; grid on;
    plot(t, m_true_B(i,:),    'k',  'LineWidth', 1.5);
    plot(t, m_lin_B(i,:),     'c--','LineWidth', 1.2);
    plot(t, m_dyn_B_aug(i,:), 'm',  'LineWidth', 1.2);
    plot(t, m_sta_B_aug(i,:), 'b--','LineWidth', 1.2);
    xline(t(k_step),'k--','Step in F_3');
    xlabel('Time [s]');
    ylabel(sprintf('m_%d [g]',i));
    title(sprintf('Scenario B – State m_%d',i));
    if i == 1
        legend('True (SDE)','Linear model', ...
               'Dyn Aug KF','Static Aug KF','Step in F_3','Location','Best');
    end
end

%% ===============================================================
% Local helper functions
% ===============================================================
function [x_lin, y_lin] = simulate_linear_ss(Ad,Bd,Cd,du,ys,N)
    % Simulate discrete linear SS model in deviation form:
    % x_{k+1} = Ad x_k + Bd du_k
    % dy_k    = Cd x_k
    % y_k     = dy_k + ys
    [nx,~] = size(Ad);
    ny     = size(Cd,1);

    x_lin = zeros(nx,N);   % deviation states
    y_lin = zeros(ny,N);   % absolute outputs

    for k = 1:N
        if k == 1
            x_prev = zeros(nx,1);
        else
            x_prev = x_lin(:,k-1);
        end
        x_lin(:,k) = Ad*x_prev + Bd*du(:,k);
        y_lin(:,k) = Cd*x_lin(:,k) + ys;
    end
end

function res = runKF_nonaug(Ad,Bd,Cd,Qx,R,dy_meas,du,ys)
    % Dynamic & static Kalman filters for non-augmented model
    % States are deviation masses x = m - m_s.
    [nx,~] = size(Ad);
    ny     = size(Cd,1);
    N      = size(dy_meas,2);

    % ---------- Dynamic KF ----------
    x_dyn = zeros(nx,N);   % deviation states
    y_dyn = zeros(ny,N);   % absolute outputs
    P_dyn = eye(nx);       % initial covariance

    for k = 1:N
        if k == 1
            x_prev = zeros(nx,1);
        else
            x_prev = x_dyn(:,k-1);
        end

        % Prediction
        x_pred = Ad*x_prev + Bd*du(:,k);
        P_pred = Ad*P_dyn*Ad' + Qx;   % Qx already contains GQG'

        % Innovation
        e  = dy_meas(:,k) - Cd*x_pred;
        Re = Cd*P_pred*Cd' + R;
        K  = P_pred*Cd'/Re;

        % Update
        x_hat = x_pred + K*e;
        P_dyn = P_pred - K*Re*K';

        x_dyn(:,k) = x_hat;
        y_dyn(:,k) = Cd*x_hat + ys;   % back to absolute output
    end

    % ---------- Static KF (steady-state gain) ----------
    [Pss,~,~] = dare(Ad',Cd',Qx,R);
    Re_ss = Cd*Pss*Cd' + R;
    Kss   = Pss*Cd'/Re_ss;

    x_sta = zeros(nx,N);
    y_sta = zeros(ny,N);

    for k = 1:N
        if k == 1
            x_prev = zeros(nx,1);
        else
            x_prev = x_sta(:,k-1);
        end
        x_pred = Ad*x_prev + Bd*du(:,k);
        e      = dy_meas(:,k) - Cd*x_pred;
        x_hat  = x_pred + Kss*e;

        x_sta(:,k) = x_hat;
        y_sta(:,k) = Cd*x_hat + ys;
    end

    res.x_dyn = x_dyn;
    res.y_dyn = y_dyn;
    res.x_sta = x_sta;
    res.y_sta = y_sta;
end

function res = runKF_aug(Ad,Bd,Cd,Qx,Qb,R,dy_meas,du,ys)
    % Dynamic & static augmented KF: state [x; b]
    % x = mass deviations, b = output bias
    [nx,nu] = size(Bd);
    ny      = size(Cd,1);
    N       = size(dy_meas,2);
    nb      = ny;

    % Augmented model
    Aa = [Ad zeros(nx,nb);
          zeros(nb,nx) eye(nb)];
    Ba = [Bd;
          zeros(nb,nu)];
    Ca = [Cd eye(ny)];

    Qa = blkdiag(Qx,Qb);

    % ---------- Dynamic augmented KF ----------
    xa_dyn = zeros(nx+nb,N);  % [x; b]
    y_dyn  = zeros(ny,N);     % absolute outputs
    b_dyn  = zeros(ny,N);     % bias estimates
    Pa     = eye(nx+nb)*10;

    for k = 1:N
        if k == 1
            xa_prev = zeros(nx+nb,1);
        else
            xa_prev = xa_dyn(:,k-1);
        end

        % Prediction
        xa_pred = Aa*xa_prev + Ba*du(:,k);
        P_pred  = Aa*Pa*Aa' + Qa;

        % Innovation
        e  = dy_meas(:,k) - Ca*xa_pred;
        Re = Ca*P_pred*Ca' + R;
        K  = P_pred*Ca'/Re;

        % Update
        xa_hat = xa_pred + K*e;
        Pa     = P_pred - K*Re*K';

        xa_dyn(:,k) = xa_hat;
        y_dyn(:,k)  = Ca*xa_hat + ys;   % absolute output
        b_dyn(:,k)  = xa_hat(nx+1:end); % bias part
    end

    % ---------- Static augmented KF ----------
    [Pss_a,~,~] = dare(Aa',Ca',Qa,R);
    Re_ss_a = Ca*Pss_a*Ca' + R;
    Kss_a   = Pss_a*Ca'/Re_ss_a;

    xa_sta = zeros(nx+nb,N);
    y_sta  = zeros(ny,N);
    b_sta  = zeros(ny,N);

    for k = 1:N
        if k == 1
            xa_prev = zeros(nx+nb,1);
        else
            xa_prev = xa_sta(:,k-1);
        end

        xa_pred = Aa*xa_prev + Ba*du(:,k);
        e       = dy_meas(:,k) - Ca*xa_pred;
        xa_hat  = xa_pred + Kss_a*e;

        xa_sta(:,k) = xa_hat;
        y_sta(:,k)  = Ca*xa_hat + ys;
        b_sta(:,k)  = xa_hat(nx+1:end);
    end

    % Pack results (only x part returned for states)
    res.x_dyn = xa_dyn(1:nx,:);
    res.x_sta = xa_sta(1:nx,:);
    res.y_dyn = y_dyn;
    res.y_sta = y_sta;
    res.b_dyn = b_dyn;
    res.b_sta = b_sta;
end
