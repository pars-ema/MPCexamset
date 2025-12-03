%% ========================================================================
% PROBLEM 8 – Offset-Free Unconstrained MPC (Azam-style)
% Using SAME linear model + noise params as in PROBLEM 6.
%
% Disturbance model (as in Problem 6 Task 3):
%   x_{k+1} = A x_k + B u_dev_k + G d_k
%   d_{k+1} = d_k                (random walk)
%   y_k     = C x_k
%
% Augmented state for KF + MPC: x_e = [x; d]
%
% Experiments:
%   A) Reference steps, true disturbances d = 0
%   B) Constant reference, step in disturbances F3,F4
%
% Requires on path:
%   - FourTankSystemModified.m
%   - FourTankSystemSensor.m
% ========================================================================

clear; clc; close all;

%% ========================================================================
% 0) Physical parameters + steady-state (same as Problem 6)
% ========================================================================

Ts = 10;                % Sampling time [s] – SAME as in Problem 6

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

fprintf("Steady-state xs (masses) from Problem 6:\n");
disp(xs);

% Steady-state heights (all 4 tanks)
h_s_abs = FourTankSystemSensor(xs,p);   % [4x1] heights in cm
output_index = [1 2];                   % we measure tanks 1 and 2
zs = h_s_abs(output_index);             % steady-state measured heights [2x1]

fprintf("Steady-state zs = [h1_s; h2_s]:\n");
disp(zs);

%% ========================================================================
% 1) Continuous linearization & discretization (same structure as Problem 6)
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
% 2) Noise covariances – EXACTLY as in Problem 6
% ========================================================================

% Process noise on disturbance channels w_k (2D: affects tanks 3 & 4)
Qw = diag([25 25]);         % same as P6

% Equivalent process noise on states: Qx = G Qw G'
Qx = G * Qw * G';           % same as P6

% Disturbance random-walk covariance for augmented model
Qd_RW = 1 * eye(nd);        % same as P6 (Qd_RW)

% Measurement noise covariance R (heights h1,h2)
R_meas = diag([4 4]);       % same as P6

%% ========================================================================
% 3) Augmented model for OFFSET-FREE MPC (x_e = [x; d])
%
%   x_{k+1} = A x_k + B u_dev_k + E d_k
%   d_{k+1} = d_k
%   y_k     = C x_k
%
% Augmented:
%   x_e(k)     = [x(k); d(k)]
%   x_e(k+1)   = A_e x_e(k) + B_e u_dev(k)
%   y(k)       = C_e x_e(k)
% ========================================================================

A_e = [A, E;
       zeros(nd,nx), eye(nd)];
B_e = [B;
       zeros(nd,nu)];
C_e = [C, zeros(ny,nd)];

nx_e = nx + nd;

%% ========================================================================
% 4) STATIC KALMAN FILTER for augmented model (as in Azam / Problem 6)
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

% Weights
Wz  = diag([1 1]);                  % output tracking weight
Wu  = 1e-4*eye(nu);                 % input magnitude weight

Qbar = kron(eye(N),Wz);             % (N*ny x N*ny)
Rbar = kron(eye(N),Wu);             % (N*nu x N*nu)

% Hessian for unconstrained QP: J = U'HU + 2*g'U
H = 2*(Gamma'*Qbar*Gamma + Rbar);

%% ========================================================================
% 6) Simulation settings (time, experiments)
% ========================================================================

Tsim = 250;                         % number of MPC steps
t    = (0:Tsim-1)*Ts;

%% ========================================================================
% 7) EXPERIMENT A: Reference steps, disturbances = 0
% ========================================================================

disp('=== EXPERIMENT A: Reference steps, disturbances = 0 ===');

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

simA = simulateOffsetFreeMPC( ...
           A,B,E,C, ...
           A_e,B_e,C_e, ...
           K_e, ...
           H,Phi,Gamma,Qbar, ...
           u_s,zs, ...
           z_ref_dev_A, ...
           d_dev_A, ...
           R_meas, ...
           Tsim,Ts);

%% ========================================================================
% 8) EXPERIMENT B: Constant reference, step in disturbances
% ========================================================================

disp('=== EXPERIMENT B: Constant reference, step in disturbances ===');

z_ref_abs_B = repmat(zs,1,Tsim);
z_ref_dev_B = zeros(2,Tsim);        % exactly at steady-state

d_dev_B = zeros(nd,Tsim);
k_step  = round(1200/Ts);           % disturbance step around 1200 s
if k_step < 1, k_step = 1; end

d_step_val = [30; -20];             % [ΔF3; ΔF4] in cm^3/s
d_dev_B(:,k_step:end) = repmat(d_step_val,1,Tsim-k_step+1);

simB = simulateOffsetFreeMPC( ...
           A,B,E,C, ...
           A_e,B_e,C_e, ...
           K_e, ...
           H,Phi,Gamma,Qbar, ...
           u_s,zs, ...
           z_ref_dev_B, ...
           d_dev_B, ...
           R_meas, ...
           Tsim,Ts);

%% ========================================================================
% 9) PLOTS
% ========================================================================

% ---------- Experiment A ----------
figure('Name','ExpA_outputs');
sgtitle('Experiment A – Reference steps, d = 0');

subplot(3,1,1); hold on; grid on;
plot(t, z_ref_abs_A(1,:), 'k--','LineWidth',1.2, 'DisplayName','h_1 ref');
plot(t, simA.y_abs(1,:),  'b',  'LineWidth',1.4, 'DisplayName','h_1');
ylabel('h_1 [cm]');
legend('Location','best');

subplot(3,1,2); hold on; grid on;
plot(t, z_ref_abs_A(2,:), 'k--','LineWidth',1.2, 'DisplayName','h_2 ref');
plot(t, simA.y_abs(2,:),  'b',  'LineWidth',1.4, 'DisplayName','h_2');
ylabel('h_2 [cm]');
legend('Location','best');

subplot(3,1,3); hold on; grid on;
plot(t, simA.u(1,:), 'b','LineWidth',1.4,'DisplayName','u_1');
plot(t, simA.u(2,:), 'r','LineWidth',1.4,'DisplayName','u_2');
yline(u_s(1),'k--','LineWidth',1.0);
yline(u_s(2),'k--','LineWidth',1.0);
ylabel('u [cm^3/s]'); xlabel('Time [s]');
legend('Location','best');

% Disturbance estimates in A (should be ~0)
figure('Name','ExpA_disturbances');
sgtitle('Experiment A – Disturbance estimates (should be ~0)');

subplot(2,1,1); hold on; grid on;
plot(t, simA.d_hat(1,:), 'b','LineWidth',1.4);
ylabel('\hat d_1 = \Delta F_3 [cm^3/s]');

subplot(2,1,2); hold on; grid on;
plot(t, simA.d_hat(2,:), 'r','LineWidth',1.4);
ylabel('\hat d_2 = \Delta F_4 [cm^3/s]');
xlabel('Time [s]');

% ---------- Experiment B ----------
figure('Name','ExpB_outputs');
sgtitle('Experiment B – Disturbance step, offset-free MPC');

subplot(3,1,1); hold on; grid on;
plot(t, z_ref_abs_B(1,:), 'k--','LineWidth',1.2,'DisplayName','h_1 ref');
plot(t, simB.y_abs(1,:),  'b',  'LineWidth',1.4,'DisplayName','h_1');
xline(t(k_step),'r--','Step in d','LabelVerticalAlignment','bottom');
ylabel('h_1 [cm]');
legend('Location','best');

subplot(3,1,2); hold on; grid on;
plot(t, z_ref_abs_B(2,:), 'k--','LineWidth',1.2,'DisplayName','h_2 ref');
plot(t, simB.y_abs(2,:),  'b',  'LineWidth',1.4,'DisplayName','h_2');
xline(t(k_step),'r--','Step in d','LabelVerticalAlignment','bottom');
ylabel('h_2 [cm]');
legend('Location','best');

subplot(3,1,3); hold on; grid on;
plot(t, simB.u(1,:), 'b','LineWidth',1.4,'DisplayName','u_1');
plot(t, simB.u(2,:), 'r','LineWidth',1.4,'DisplayName','u_2');
xline(t(k_step),'r--','Step in d');
yline(u_s(1),'k--','LineWidth',1.0);
yline(u_s(2),'k--','LineWidth',1.0);
ylabel('u [cm^3/s]'); xlabel('Time [s]');
legend('Location','best');

% Disturbances true vs estimated (Experiment B)
figure('Name','ExpB_disturbances');
sgtitle('Experiment B – True vs estimated disturbances');

subplot(2,1,1); hold on; grid on;
plot(t, d_dev_B(1,:),    'k--','LineWidth',1.4,'DisplayName','True \DeltaF_3');
plot(t, simB.d_hat(1,:), 'b',  'LineWidth',1.4,'DisplayName','\hat{\DeltaF}_3');
xline(t(k_step),'r--','Step');
ylabel('\DeltaF_3 [cm^3/s]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(t, d_dev_B(2,:),    'k--','LineWidth',1.4,'DisplayName','True \DeltaF_4');
plot(t, simB.d_hat(2,:), 'b',  'LineWidth',1.4,'DisplayName','\hat{\DeltaF}_4');
xline(t(k_step),'r--','Step');
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
function sim = simulateOffsetFreeMPC( ...
            A,B,E,C, ...
            A_e,B_e,C_e, ...
            K_e, ...
            H,Phi,Gamma,Qbar, ...
            u_s,zs, ...
            z_ref_dev, ...
            d_dev, ...
            R_meas, ...
            Tsim,Ts)
%SIMULATEOFFSETFREEMPC
% Closed-loop MPC + static augmented KF.
%
% A,B,E,C      : deviation model matrices
% A_e,B_e,C_e  : augmented model (x;d)
% K_e          : static KF gain (from DARE)
% H,Phi,Gamma,Qbar: MPC matrices for augmented model
% u_s,zs       : steady-state input and output
% z_ref_dev    : deviation reference [2 x Tsim]
% d_dev        : true disturbance deviations [2 x Tsim]
% R_meas       : measurement noise covariance
% Tsim,Ts      : number of steps and sampling time

    [nx,nu] = size(B);
    nd      = size(E,2);
    ny      = size(C,1);

    nx_e    = nx + nd;
    N       = size(Phi,1)/ny;

    % Allocate
    x_dev   = zeros(nx, Tsim+1);    % true state deviation
    d_true  = zeros(nd, Tsim+1);    % true disturbance dev
    y_true  = zeros(ny, Tsim);      % true output dev
    y_abs   = zeros(ny, Tsim);      % absolute heights
    y_meas  = zeros(ny, Tsim);      % measured dev

    x_hat_e = zeros(nx_e, Tsim+1);  % KF augmented state [x; d]
    u_dev   = zeros(nu, Tsim);      % control dev
    u_abs   = zeros(nu, Tsim);      % absolute inputs

    d_hat   = zeros(nd, Tsim);      % estimated disturbance dev

    % Initial conditions
    x_dev(:,1)   = zeros(nx,1);     % start at steady state
    d_true(:,1)  = d_dev(:,1);      % initial disturbance
    x_hat_e(:,1) = zeros(nx_e,1);   % initial KF estimate [0;0]

    sigma_y = sqrt(diag(R_meas));   % measurement noise std dev

    % Closed-loop simulation
    for k = 1:Tsim

        % --- True plant update (for k>1) ---
        if k > 1
            x_dev(:,k) = A*x_dev(:,k-1) + B*u_dev(:,k-1) + E*d_true(:,k-1);
        end
        y_true(:,k) = C*x_dev(:,k);
        y_abs(:,k)  = zs + y_true(:,k);

        % Measurement with noise
        noise_k = sigma_y .* randn(ny,1);
        y_meas(:,k) = y_true(:,k) + noise_k;

        % True disturbance dev
        if k < Tsim
            d_true(:,k+1) = d_dev(:,k+1);
        end

        % --- Static augmented Kalman filter ---
        innov        = y_meas(:,k) - C_e*x_hat_e(:,k);
        x_hat_e(:,k) = x_hat_e(:,k) + K_e*innov;

        x_hat_k = x_hat_e(1:nx,k);
        d_hat_k = x_hat_e(nx+1:end,k);
        d_hat(:,k) = d_hat_k;   %#ok<NASGU>  (logged below anyway)

        % --- MPC: build stacked reference over horizon ---
        zref_stack = zeros(ny*N,1);
        for j = 0:N-1
            idx = min(k+j, Tsim);
            zref_stack(j*ny+(1:ny)) = z_ref_dev(:,idx);
        end

        % Quadratic cost gradient
        g = 2*Gamma'*Qbar*(Phi*x_hat_e(:,k) - zref_stack);

        % Unconstrained optimal moves
        U_seq = -H \ g;
        u_dev(:,k) = U_seq(1:nu);
        u_abs(:,k) = u_s + u_dev(:,k);

        % KF time update
        if k < Tsim
            x_hat_e(:,k+1) = A_e*x_hat_e(:,k) + B_e*u_dev(:,k);
        end
    end

    % Output struct
    sim = struct();
    sim.x_dev = x_dev;
    sim.y_dev = y_true;
    sim.y_abs = y_abs;
    sim.y_meas = y_meas;
    sim.u_dev = u_dev;
    sim.u     = u_abs;
    sim.d_true = d_true(:,1:Tsim);
    sim.d_hat  = x_hat_e(nx+1:end,1:Tsim);
end
