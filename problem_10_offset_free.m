%% ========================================================================
% PROBLEM 10 – Offset-Free MPC on LINEAR vs NONLINEAR Models
%
% Controller:
%   - SAME constrained offset-free MPC as in Problem 9
%   - SAME augmented static Kalman Filter (x + d) as in Problem 6/9
%
% Plants compared:
%   - Linear plant:      discrete-time linear model from linearization
%   - Nonlinear plant:   FourTankSystemModified + FourTankSystemSensor
%
% Experiments:
%   A) Reference steps, disturbances d_dev = 0
%   B) SAME reference steps as A, PLUS disturbance step later
%
% QP solved with problem_7.m (must be on MATLAB path):
%   [U,info] = problem_7(H,g,l,u,A,bl,bu,xinit)
%
% Requires on path:
%   - FourTankSystemModified.m
%   - FourTankSystemSensor.m
%   - problem_7.m
% ========================================================================

clear; clc; close all;

%% ========================================================================
% 0) Physical parameters + steady-state (same as Problem 9, Ts = 10 s)
% ========================================================================

Ts = 10;                % Sampling time [s]

a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;        % [cm^2]
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327;% [cm^2]
g  = 981;                                                   % [cm/s^2]
rho = 1;                                                    % [g/cm^3]
gamma1 = 0.58;
gamma2 = 0.72;

p = [a1;a2;a3;a4;A1;A2;A3;A4;g;gamma1;gamma2;rho];

% Operating point
u_s = [300;300];          % [cm^3/s] (F1, F2)
d_s = [250;250];          % [cm^3/s] (F3, F4) – only used for nonlinear model

% Steady-state masses xs (for nonlinear model)
xs_guess = 500*ones(4,1);
opts = optimoptions('fsolve','Display','none');
xs = fsolve(@(x)FourTankSystemModified(0,x,u_s,d_s,p), xs_guess, opts);

fprintf("Steady-state xs (masses):\n");
disp(xs);

% Steady-state heights (all 4 tanks)
h_s_abs = FourTankSystemSensor(xs,p);   % [4x1] heights in cm
output_index = [1 2];                   % we measure tanks 1 and 2
zs = h_s_abs(output_index);             % steady-state measured heights [2x1]

fprintf("Steady-state zs = [h1_s; h2_s]:\n");
disp(zs);

%% ========================================================================
% 1) Continuous linearization & discretization (same structure as P9)
% ========================================================================

nx = 4; nu = 2; nd = 2; ny = 2;

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

% Discretize [A_c, B_c, E_c] with ZOH
M = [A_c, [B_c E_c];
     zeros(nu+nd, nx+nu+nd)];
Md = expm(M*Ts);

Ad = Md(1:nx, 1:nx);
BEd = Md(1:nx, nx+1:end);
Bd  = BEd(:, 1:nu);
Ed  = BEd(:, nu+1:end);     % [nx x nd]

Cd = C_c;

A = Ad;
B = Bd;
E = Ed;
C = Cd;

G = E;                      % disturbance / noise input matrix

%% ========================================================================
% 2) Noise covariances – SAME as Problem 6/9
% ========================================================================

% Process noise on disturbance channels w_k (2D: affects tanks 3 & 4)
Qw = diag([25 25]);         % same as P6/P9

% Equivalent process noise on states: Qx = G Qw G'
Qx = G * Qw * G';           % same as P6/P9

% Disturbance random-walk covariance for augmented model
Qd_RW = 1 * eye(nd);        % same as P6/P9

% Measurement noise covariance R (heights h1,h2)
R_meas = diag([4 4]);       % std 2 cm on each level

%% ========================================================================
% 3) Augmented model for OFFSET-FREE MPC (x_e = [x; d])
% ========================================================================

A_e = [A, E;
       zeros(nd,nx), eye(nd)];
B_e = [B;
       zeros(nd,nu)];
C_e = [C, zeros(ny,nd)];

nx_e = nx + nd;

%% ========================================================================
% 4) STATIC KALMAN FILTER for augmented model (as in Azam / P6/P9)
% ========================================================================

Q_e = blkdiag(Qx, Qd_RW);           % process covariance for [x; d]

[P_e,~,~] = dare(A_e',C_e',Q_e,R_meas);
P_e = 0.5*(P_e+P_e');               % enforce symmetry

Re = C_e*P_e*C_e' + R_meas;
K_e = P_e * C_e' / Re;              % [nx_e x ny] Kalman gain

fprintf("Static augmented KF gain norm: %.4f\n", norm(K_e));

%% ========================================================================
% 5) MPC prediction matrices (for augmented model)
% ========================================================================

N = 15;                             % prediction horizon

[Phi_x_e, Gamma_x_e] = buildPredictionMatrices(A_e,B_e,N);
Phi   = kron(eye(N),C_e) * Phi_x_e;      % (N*ny x nx_e)
Gamma = kron(eye(N),C_e) * Gamma_x_e;   % (N*ny x N*nu)

% Weights (same as Problem 9)
Wz  = diag([1 1]);                  % output tracking weight
Wu  = 1e-4*eye(nu);                 % input magnitude weight

Qbar = kron(eye(N),Wz);             % (N*ny x N*ny)
Rbar = kron(eye(N),Wu);             % (N*nu x N*nu)

% Hessian for constrained QP: 0.5*U'HU + f'U
H = 2*(Gamma'*Qbar*Gamma + Rbar);

%% ========================================================================
% 6) Input amplitude constraints (absolute pumps)
%     0 <= u <= 600  (cm^3/s)
% implemented via bounds on U_dev = [u_dev(k)...]
% ========================================================================

u_min_abs = [0; 0];
u_max_abs = [600; 600];

% Deviation constraints: u_dev = u - u_s
lb_u_dev = u_min_abs - u_s;
ub_u_dev = u_max_abs - u_s;

% Bounds on stacked U_dev = [u_dev(k); ...; u_dev(k+N-1)]
lb = repmat(lb_u_dev, N, 1);    % (N*nu x 1)
ub = repmat(ub_u_dev, N, 1);    % (N*nu x 1)

%% ========================================================================
% 7) Simulation settings (time, experiments)
% ========================================================================

Tsim = 300;                         % number of MPC steps
t    = (0:Tsim-1)*Ts;

rng(0);                             % reproducibility for measurement noise

%% ========================================================================
% 8) EXPERIMENT A: Reference steps, disturbances = 0
% ========================================================================

disp('=== PROBLEM 10 – EXPERIMENT A: Reference steps, d_dev = 0 ===');

t_change_h1 = [0 600 1600];
values_h1   = [zs(1) zs(1)*1.3 zs(1)*0.8];

t_change_h2 = [0 600 1600];
values_h2   = [zs(2) zs(2)*1.3 zs(2)*1.0];

[h1_ref_A, h2_ref_A] = generate_stair_references( ...
    t, t_change_h1, values_h1, ...
       t_change_h2, values_h2);

z_ref_abs_A = [h1_ref_A; h2_ref_A];
z_ref_dev_A = z_ref_abs_A - zs;     % deviation reference

% Disturbances dev in Experiment A: zero
d_dev_A = zeros(nd,Tsim);

% --- A1: Linear plant (reuse linear model internally) ---
simA_lin = simulateOffsetFreeMPC_constrained_linear( ...
           A,B,E,C, ...
           A_e,B_e,C_e, ...
           K_e, ...
           H,Phi,Gamma,Qbar, ...
           u_s,zs, ...
           z_ref_dev_A, ...
           d_dev_A, ...
           R_meas, ...
           lb,ub, ...
           Tsim,Ts);

% --- A2: Nonlinear plant (FourTankSystemModified + Sensor) ---
simA_nl = simulateOffsetFreeMPC_constrained_nonlinear( ...
           A,B,E,C, ...
           A_e,B_e,C_e, ...
           K_e, ...
           H,Phi,Gamma,Qbar, ...
           u_s,zs, ...
           z_ref_dev_A, ...
           d_dev_A, ...
           R_meas, ...
           lb,ub, ...
           Tsim,Ts, ...
           p,d_s,xs,output_index);

%% ========================================================================
% 9) EXPERIMENT B: SAME reference steps as A, PLUS disturbance step later
% ========================================================================

disp('=== PROBLEM 10 – EXPERIMENT B: Reference steps + disturbance step ===');

% Same reference steps as in A
z_ref_abs_B = z_ref_abs_A;
z_ref_dev_B = z_ref_dev_A;

% Disturbance step (later)
d_dev_B = zeros(nd,Tsim);
k_step = round(2000/Ts);   % when to inject disturbance step

d_step_val = [30; -20];    % [ΔF3; ΔF4]
if k_step < 1, k_step = 1; end
d_dev_B(:,k_step:end) = repmat(d_step_val,1,Tsim-k_step+1);

% --- B1: Linear plant ---
simB_lin = simulateOffsetFreeMPC_constrained_linear( ...
           A,B,E,C, ...
           A_e,B_e,C_e, ...
           K_e, ...
           H,Phi,Gamma,Qbar, ...
           u_s,zs, ...
           z_ref_dev_B, ...
           d_dev_B, ...
           R_meas, ...
           lb,ub, ...
           Tsim,Ts);

% --- B2: Nonlinear plant ---
simB_nl = simulateOffsetFreeMPC_constrained_nonlinear( ...
           A,B,E,C, ...
           A_e,B_e,C_e, ...
           K_e, ...
           H,Phi,Gamma,Qbar, ...
           u_s,zs, ...
           z_ref_dev_B, ...
           d_dev_B, ...
           R_meas, ...
           lb,ub, ...
           Tsim,Ts, ...
           p,d_s,xs,output_index);

%% ========================================================================
% 10) PLOTS – Linear vs Nonlinear (Problem 10)
% ========================================================================

% ---------- Experiment A ----------
figure('Name','P10_ExpA_outputs');
sgtitle('Problem 10 – Exp A: Reference steps, linear vs nonlinear plant');

subplot(3,1,1); hold on; grid on;
plot(t, z_ref_abs_A(1,:),      'k--','LineWidth',1.2, 'DisplayName','h_1 ref');
plot(t, simA_lin.y_abs(1,:),   'b',  'LineWidth',1.4, 'DisplayName','h_1 linear');
plot(t, simA_nl.y_abs(1,:),    'r',  'LineWidth',1.4, 'DisplayName','h_1 nonlinear');
ylabel('h_1 [cm]');
legend('Location','best');

subplot(3,1,2); hold on; grid on;
plot(t, z_ref_abs_A(2,:),      'k--','LineWidth',1.2, 'DisplayName','h_2 ref');
plot(t, simA_lin.y_abs(2,:),   'b',  'LineWidth',1.4, 'DisplayName','h_2 linear');
plot(t, simA_nl.y_abs(2,:),    'r',  'LineWidth',1.4, 'DisplayName','h_2 nonlinear');
ylabel('h_2 [cm]');
legend('Location','best');

subplot(3,1,3); hold on; grid on;
plot(t, simA_lin.u(1,:), 'b','LineWidth',1.4,'DisplayName','u_1 linear');
plot(t, simA_lin.u(2,:), 'b:','LineWidth',1.2,'DisplayName','u_2 linear');
plot(t, simA_nl.u(1,:),  'r','LineWidth',1.4,'DisplayName','u_1 nonlinear');
plot(t, simA_nl.u(2,:),  'r:','LineWidth',1.2,'DisplayName','u_2 nonlinear');
yline(u_s(1),'k--','LineWidth',1.0);
yline(u_s(2),'k--','LineWidth',1.0);
yline(0,'k:');
yline(600,'k:');
ylabel('u [cm^3/s]'); xlabel('Time [s]');
legend('Location','best');

% Disturbance estimates in A (should be ~0)
figure('Name','P10_ExpA_disturbances');
sgtitle('Problem 10 – Exp A: Disturbance estimates (linear vs nonlinear)');

subplot(2,1,1); hold on; grid on;
plot(t, simA_lin.d_hat(1,:), 'b','LineWidth',1.4,'DisplayName','\hat d_1 linear');
plot(t, simA_nl.d_hat(1,:),  'r','LineWidth',1.4,'DisplayName','\hat d_1 nonlinear');
ylabel('\hat d_1 = \Delta F_3 [cm^3/s]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(t, simA_lin.d_hat(2,:), 'b','LineWidth',1.4,'DisplayName','\hat d_2 linear');
plot(t, simA_nl.d_hat(2,:),  'r','LineWidth',1.4,'DisplayName','\hat d_2 nonlinear');
ylabel('\hat d_2 = \Delta F_4 [cm^3/s]');
xlabel('Time [s]');
legend('Location','best');

% ---------- Experiment B ----------
figure('Name','P10_ExpB_outputs');
sgtitle('Problem 10 – Exp B: Disturbance step, linear vs nonlinear plant');

subplot(3,1,1); hold on; grid on;
plot(t, z_ref_abs_B(1,:),      'k--','LineWidth',1.2,'DisplayName','h_1 ref');
plot(t, simB_lin.y_abs(1,:),   'b',  'LineWidth',1.4,'DisplayName','h_1 linear');
plot(t, simB_nl.y_abs(1,:),    'r',  'LineWidth',1.4,'DisplayName','h_1 nonlinear');
xline(t(k_step),'m--','Step in d','LabelVerticalAlignment','bottom');
ylabel('h_1 [cm]');
legend('Location','best');

subplot(3,1,2); hold on; grid on;
plot(t, z_ref_abs_B(2,:),      'k--','LineWidth',1.2,'DisplayName','h_2 ref');
plot(t, simB_lin.y_abs(2,:),   'b',  'LineWidth',1.4,'DisplayName','h_2 linear');
plot(t, simB_nl.y_abs(2,:),    'r',  'LineWidth',1.4,'DisplayName','h_2 nonlinear');
xline(t(k_step),'m--','Step in d','LabelVerticalAlignment','bottom');
ylabel('h_2 [cm]');
legend('Location','best');

subplot(3,1,3); hold on; grid on;
plot(t, simB_lin.u(1,:), 'b','LineWidth',1.4,'DisplayName','u_1 linear');
plot(t, simB_lin.u(2,:), 'b:','LineWidth',1.2,'DisplayName','u_2 linear');
plot(t, simB_nl.u(1,:),  'r','LineWidth',1.4,'DisplayName','u_1 nonlinear');
plot(t, simB_nl.u(2,:),  'r:','LineWidth',1.2,'DisplayName','u_2 nonlinear');
xline(t(k_step),'m--','Step in d');
yline(u_s(1),'k--','LineWidth',1.0);
yline(u_s(2),'k--','LineWidth',1.0);
yline(0,'k:');
yline(600,'k:');
ylabel('u [cm^3/s]'); xlabel('Time [s]');
legend('Location','best');

% Disturbances true vs estimated (Experiment B)
figure('Name','P10_ExpB_disturbances');
sgtitle('Problem 10 – Exp B: True vs estimated disturbances');

subplot(2,1,1); hold on; grid on;
plot(t, d_dev_B(1,:),        'k--','LineWidth',1.4,'DisplayName','True \DeltaF_3');
plot(t, simB_lin.d_hat(1,:), 'b',  'LineWidth',1.4,'DisplayName','\hat{\DeltaF}_3 linear');
plot(t, simB_nl.d_hat(1,:),  'r',  'LineWidth',1.4,'DisplayName','\hat{\DeltaF}_3 nonlinear');
xline(t(k_step),'m--','Step');
ylabel('\DeltaF_3 [cm^3/s]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(t, d_dev_B(2,:),        'k--','LineWidth',1.4,'DisplayName','True \DeltaF_4');
plot(t, simB_lin.d_hat(2,:), 'b',  'LineWidth',1.4,'DisplayName','\hat{\DeltaF}_4 linear');
plot(t, simB_nl.d_hat(2,:),  'r',  'LineWidth',1.4,'DisplayName','\hat{\DeltaF}_4 nonlinear');
xline(t(k_step),'m--','Step');
ylabel('\DeltaF_4 [cm^3/s]'); xlabel('Time [s]');
legend('Location','best');

%% ========================================================================
% ===================== HELPER FUNCTIONS ==================================
% ========================================================================

function [Phi_x, Gamma_x] = buildPredictionMatrices(A,B,N)
%BUILD PREDICTION MATRICES FOR:
%   x_{k+1} = A x_k + B u_k
%
% Phi_x   : [N*nx x nx]
% Gamma_x : [N*nx x N*nu]

    [nx,nu] = size(B);
    Phi_x   = zeros(N*nx, nx);
    Gamma_x = zeros(N*nx, N*nu);

    A_power = eye(nx);

    for i = 1:N
        A_power = A * A_power;          % A^i
        rows = (i-1)*nx + (1:nx);
        Phi_x(rows,:) = A_power;

        % input influence
        A_j = eye(nx);
        for j = 1:i
            cols = (j-1)*nu + (1:nu);
            Gamma_x(rows,cols) = A_j * B;
            A_j = A * A_j;
        end
    end
end

% -------------------------------------------------------------------------
function sim = simulateOffsetFreeMPC_constrained_linear( ...
            A,B,E,C, ...
            A_e,B_e,C_e, ...
            K_e, ...
            H,Phi,Gamma,Qbar, ...
            u_s,zs, ...
            z_ref_dev, ...
            d_dev, ...
            R_meas, ...
            lb,ub, ...
            Tsim,Ts)
% Closed-loop constrained MPC + static augmented KF, LINEAR plant.

    [nx,nu] = size(B);
    nd      = size(E,2);
    ny      = size(C,1);

    nx_e    = nx + nd;
    N       = size(Phi,1)/ny;

    % Allocate
    x_dev   = zeros(nx, Tsim+1);    % true state deviation (linear)
    y_true  = zeros(ny, Tsim);      % true output dev
    y_abs   = zeros(ny, Tsim);      % absolute heights
    y_meas  = zeros(ny, Tsim);      % measured dev

    x_hat_e = zeros(nx_e, Tsim+1);  % KF augmented state [x; d]
    u_dev   = zeros(nu, Tsim);      % control dev
    u_abs   = zeros(nu, Tsim);      % absolute inputs

    d_hat   = zeros(nd, Tsim);      % estimated disturbance dev

    % Initial conditions
    x_dev(:,1)   = zeros(nx,1);     % start at steady state
    x_hat_e(:,1) = zeros(nx_e,1);   % initial KF estimate [0;0]

    sigma_y = sqrt(diag(R_meas));   % measurement noise std dev

    % QP static parts
    l = lb;
    u = ub;
    Aqp = [];       % no linear constraints, only bounds
    bl  = [];
    bu  = [];

    % Closed-loop simulation
    for k = 1:Tsim

        % --- True LINEAR plant update (for k>1) ---
        if k > 1
            x_dev(:,k) = A*x_dev(:,k-1) + B*u_dev(:,k-1) + E*d_dev(:,k-1);
        end
        y_true(:,k) = C*x_dev(:,k);
        y_abs(:,k)  = zs + y_true(:,k);

        % Measurement with noise
        noise_k = sigma_y .* randn(ny,1);
        y_meas(:,k) = y_true(:,k) + noise_k;

        % --- Static augmented Kalman filter ---
        innov        = y_meas(:,k) - C_e*x_hat_e(:,k);
        x_hat_e(:,k) = x_hat_e(:,k) + K_e*innov;

        % Store estimated disturbances
        d_hat_k = x_hat_e(nx+1:end,k);
        d_hat(:,k) = d_hat_k;

        % --- MPC: build stacked reference over horizon ---
        zref_stack = zeros(ny*N,1);
        for j = 0:N-1
            idx = min(k+j, Tsim);
            zref_stack(j*ny+(1:ny)) = z_ref_dev(:,idx);
        end

        % Quadratic cost f vector (for 0.5*U'HU + f'U)
        f = 2*Gamma'*Qbar*(Phi*x_hat_e(:,k) - zref_stack);

        % Initial guess
        xinit = zeros(N*nu,1);

        % Solve QP with bounds only
        [U_seq, info] = problem_7(H, f, l, u, Aqp, bl, bu, xinit); %#ok<NASGU>

        % First control move
        u_dev(:,k) = U_seq(1:nu);
        u_abs(:,k) = u_s + u_dev(:,k);

        % KF time update
        if k < Tsim
            x_hat_e(:,k+1) = A_e*x_hat_e(:,k) + B_e*u_dev(:,k);
        end
    end

    % Output struct
    sim = struct();
    sim.x_dev  = x_dev;
    sim.y_dev  = y_true;
    sim.y_abs  = y_abs;
    sim.y_meas = y_meas;
    sim.u_dev  = u_dev;
    sim.u      = u_abs;
    sim.d_true = d_dev;
    sim.d_hat  = x_hat_e(nx+1:end,1:Tsim);
end

% -------------------------------------------------------------------------
function sim = simulateOffsetFreeMPC_constrained_nonlinear( ...
            A,B,E,C, ...
            A_e,B_e,C_e, ...
            K_e, ...
            H,Phi,Gamma,Qbar, ...
            u_s,zs, ...
            z_ref_dev, ...
            d_dev, ...
            R_meas, ...
            lb,ub, ...
            Tsim,Ts, ...
            p,d_s,xs,output_index)
% Closed-loop constrained MPC + static augmented KF, NONLINEAR plant:
%   x_dot = FourTankSystemModified(t, x, u_abs, d_abs, p)
%   y     = FourTankSystemSensor(x,p);  (heights)
% with measurement noise and same linear-KF + linear-model-based MPC.

    [nx,nu] = size(B);
    nd      = size(E,2);
    ny      = numel(output_index);  % measured outputs (1 and 2)

    nx_e    = nx + nd;
    N       = size(Phi,1)/ny;

    % Allocate
    x_nl    = zeros(nx, Tsim+1);    % true nonlinear masses
    y_true  = zeros(ny, Tsim);      % true output dev (heights)
    y_abs   = zeros(ny, Tsim);      % absolute heights
    y_meas  = zeros(ny, Tsim);      % measured dev

    x_hat_e = zeros(nx_e, Tsim+1);  % KF augmented state [x; d] (linear model)
    u_dev   = zeros(nu, Tsim);      % control dev
    u_abs   = zeros(nu, Tsim);      % absolute inputs

    d_hat   = zeros(nd, Tsim);      % estimated disturbance dev

    % Initial conditions: start at steady state xs
    x_nl(:,1)   = xs;
    x_hat_e(:,1) = zeros(nx_e,1);   % linear KF believes we start at dev 0

    sigma_y = sqrt(diag(R_meas));   % measurement noise std dev (applied to h1,h2)

    % QP static parts
    l = lb;
    u = ub;
    Aqp = [];       % no linear constraints, only bounds
    bl  = [];
    bu  = [];

    % Closed-loop simulation
    for k = 1:Tsim

        % --- Nonlinear plant output ---
        h_all = FourTankSystemSensor(x_nl(:,k), p);  % all 4 heights
        y_abs_k = h_all(output_index);               % measured tanks 1 & 2
        y_true(:,k) = y_abs_k - zs;                  % deviation
        y_abs(:,k)  = y_abs_k;

        % Measurement with noise
        noise_k = sigma_y .* randn(ny,1);
        y_meas(:,k) = y_true(:,k) + noise_k;

        % --- Static augmented Kalman filter (linear model) ---
        innov        = y_meas(:,k) - C_e*x_hat_e(:,k);
        x_hat_e(:,k) = x_hat_e(:,k) + K_e*innov;

        % Store estimated disturbances
        d_hat_k = x_hat_e(nx+1:end,k);
        d_hat(:,k) = d_hat_k;

        % --- MPC: build stacked reference over horizon ---
        zref_stack = zeros(ny*N,1);
        for j = 0:N-1
            idx = min(k+j, Tsim);
            zref_stack(j*ny+(1:ny)) = z_ref_dev(:,idx);
        end

        % Quadratic cost f vector (for 0.5*U'HU + f'U)
        f = 2*Gamma'*Qbar*(Phi*x_hat_e(:,k) - zref_stack);

        % Initial guess
        xinit = zeros(N*nu,1);

        % Solve QP with bounds only
        [U_seq, info] = problem_7(H, f, l, u, Aqp, bl, bu, xinit); %#ok<NASGU>

        % First control move
        u_dev(:,k) = U_seq(1:nu);
        u_abs(:,k) = u_s + u_dev(:,k);

        % --- Nonlinear plant update (Euler integration over Ts) ---
        d_abs_k = d_s + d_dev(:,k);    % absolute disturbances F3,F4
        dxdt = FourTankSystemModified(0, x_nl(:,k), u_abs(:,k), d_abs_k, p);
        x_nl(:,k+1) = x_nl(:,k) + Ts * dxdt;

        % KF time update (linear model)
        if k < Tsim
            x_hat_e(:,k+1) = A_e*x_hat_e(:,k) + B_e*u_dev(:,k);
        end
    end

    % Output struct
    sim = struct();
    sim.x_nl  = x_nl;
    sim.y_dev = y_true;
    sim.y_abs = y_abs;
    sim.y_meas = y_meas;
    sim.u_dev  = u_dev;
    sim.u      = u_abs;
    sim.d_true = d_dev;
    sim.d_hat  = x_hat_e(nx+1:end,1:Tsim);
end

% -------------------------------------------------------------------------
function [h1_ref, h2_ref] = generate_stair_references( ...
                                    t, t_change_h1, values_h1, ...
                                       t_change_h2, values_h2)
% Simple staircase generator (same structure as in Problem 8/9)

    Nt = numel(t);
    h1_ref = zeros(1,Nt);
    h2_ref = zeros(1,Nt);

    for k = 1:Nt
        tk = t(k);

        % h1
        idx1 = find(tk >= t_change_h1, 1, 'last');
        if isempty(idx1), idx1 = 1; end
        idx1 = min(idx1, numel(values_h1));
        h1_ref(k) = values_h1(idx1);

        % h2
        idx2 = find(tk >= t_change_h2, 1, 'last');
        if isempty(idx2), idx2 = 1; end
        idx2 = min(idx2, numel(values_h2));
        h2_ref(k) = values_h2(idx2);
    end
end
