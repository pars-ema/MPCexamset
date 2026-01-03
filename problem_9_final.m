%% ========================================================================
% PROBLEM 9 – Input Constrained Offset-Free MPC (Azam-style)
% + overlay with NONLINEAR SDE plant (Euler–Maruyama, closed loop)
% + SIMPLE LaTeX that always works:
%       - Legends: plain text (no LaTeX)
%       - Axis labels: LaTeX only where needed (Interpreter='latex')
% + Step marker (vertical dashed red line) + movable text "step in d"
% + Save all figures as vector PDFs to figures/problem_9
%
% Requires on path:
%   - FourTankSystemModified.m
%   - FourTankSystemSensor.m
%   - ScalarStdWienerProcess.m
%   - problem_7.m
% ========================================================================

clear; clc; close all;

%% ========================================================================
% Global figure style (bigger + report-friendly)
% ========================================================================
figDir = fullfile('figures','problem_9');
if ~exist(figDir,'dir'), mkdir(figDir); end

set(groot,'defaultFigureColor','w');
set(groot,'defaultFigurePaperPositionMode','auto');
set(groot,'defaultFigurePosition',[100 100 1350 780]);   % bigger

% Keep interpreters robust:
set(groot,'defaultTextInterpreter','none');
set(groot,'defaultLegendInterpreter','none');
set(groot,'defaultAxesTickLabelInterpreter','none');

% Fonts
set(groot,'defaultAxesFontSize',13);
set(groot,'defaultLegendFontSize',12);

rng(0); % reproducibility for measurement noise

%% ========================================================================
% 0) Physical parameters + steady-state
% ========================================================================
Ts = 10;  % sampling time [s]

a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;        % [cm^2]
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327;% [cm^2]
g  = 981;                                                   % [cm/s^2]
rho = 1;                                                    % [g/cm^3]
gamma1 = 0.58;
gamma2 = 0.72;

p = [a1;a2;a3;a4;A1;A2;A3;A4;g;gamma1;gamma2;rho];

u_s = [300;300];          % [cm^3/s] (F1, F2)
d_s = [250;250];          % [cm^3/s] (F3, F4) – for nonlinear model

xs_guess = 500*ones(4,1);
opts = optimoptions('fsolve','Display','none');
xs = fsolve(@(x)FourTankSystemModified(0,x,u_s,d_s,p), xs_guess, opts);

fprintf("Steady-state xs (masses):\n"); disp(xs);

h_s_abs = FourTankSystemSensor(xs,p);   % [4x1] heights [cm]
output_index = [1 2];                   % measure tanks 1 and 2
zs = h_s_abs(output_index);             % steady-state measured heights [2x1]
fprintf("Steady-state zs = [h1_s; h2_s]:\n"); disp(zs);

%% ========================================================================
% 1) Linearization & discretization
% ========================================================================
nx = 4; nu = 2; nd = 2; ny = 2;

% unpack
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

A_c = rho * [
   -beta1,  0,       beta3,   0;
         0, -beta2,      0,  beta4;
         0,      0, -beta3,   0;
         0,      0,      0, -beta4
];

B_c = rho * [
    gamma1,        0;
         0,   gamma2;
         0, 1-gamma2;
   1-gamma1,       0
];

E_c = [
    0, 0;
    0, 0;
    1, 0;
    0, 1
];

C_c = [
    1/(rho*A1), 0, 0, 0;
    0, 1/(rho*A2), 0, 0
];

% ZOH discretization of [A_c, B_c, E_c]
M  = [A_c, [B_c E_c];
      zeros(nu+nd, nx+nu+nd)];
Md = expm(M*Ts);

Ad  = Md(1:nx, 1:nx);
BEd = Md(1:nx, nx+1:end);
Bd  = BEd(:, 1:nu);
Ed  = BEd(:, nu+1:end);

A = Ad; B = Bd; E = Ed; C = C_c;

%% ========================================================================
% 2) Noise covariances (same style as before)
% ========================================================================
Qw     = diag([25 25]);
Qx     = E * Qw * E';
Qd_RW  = 1 * eye(nd);
R_meas = diag([4 4]);

%% ========================================================================
% 3) Augmented model for offset-free MPC: x_e = [x; d]
% ========================================================================
A_e = [A, E;
       zeros(nd,nx), eye(nd)];
B_e = [B;
       zeros(nd,nu)];
C_e = [C, zeros(ny,nd)];
nx_e = nx + nd;

%% ========================================================================
% 4) Static KF for augmented model
% ========================================================================
Q_e = blkdiag(Qx, Qd_RW);
[P_e,~,~] = dare(A_e', C_e', Q_e, R_meas);
P_e = 0.5*(P_e+P_e');

Re = C_e*P_e*C_e' + R_meas;
K_e = P_e * C_e' / Re;
fprintf("Static augmented KF gain norm: %.4f\n", norm(K_e));

%% ========================================================================
% 5) MPC prediction matrices
% ========================================================================
N = 15;
[Phi_x_e, Gamma_x_e] = buildPredictionMatrices(A_e,B_e,N);
Phi   = kron(eye(N),C_e) * Phi_x_e;
Gamma = kron(eye(N),C_e) * Gamma_x_e;

Wz  = diag([1 1]);
Wu  = 1e-4*eye(nu);
Qbar = kron(eye(N),Wz);
Rbar = kron(eye(N),Wu);

% QP form: 0.5*U'HU + f'U  (so our H below is "2*(...)" and f is gradient)
H = 2*(Gamma'*Qbar*Gamma + Rbar);

%% ========================================================================
% 6) Input amplitude constraints (absolute pumps): 0 <= u <= 600
% bounds on U_dev (stacked)
% ========================================================================
u_min_abs = [0; 0];
u_max_abs = [600; 600];

lb_u_dev = u_min_abs - u_s;
ub_u_dev = u_max_abs - u_s;

lb = repmat(lb_u_dev, N, 1);
ub = repmat(ub_u_dev, N, 1);

%% ========================================================================
% 7) Simulation settings + experiments
% ========================================================================
Tsim = 300;
t    = (0:Tsim-1)*Ts;

% Experiment A references (use exactly what you wrote)
t_change_h1 = [0 600 1600];
values_h1   = [zs(1) zs(1)*1.3 zs(1)*0.8];

t_change_h2 = [0 600 1300];
values_h2   = [zs(2) zs(2)*1.3 zs(2)*1.0];

[h1_ref_A, h2_ref_A] = generate_stair_references(t, t_change_h1, values_h1, t_change_h2, values_h2);
z_ref_abs_A = [h1_ref_A; h2_ref_A];
z_ref_dev_A = z_ref_abs_A - zs;

% Experiment B uses same refs as A
z_ref_abs_B = z_ref_abs_A;
z_ref_dev_B = z_ref_dev_A;

% Disturbance sequences
d_dev_A = zeros(nd,Tsim);

d_dev_B = zeros(nd,Tsim);
k_step = round(2000/Ts);
k_step = max(1, min(Tsim, k_step));

d_step_val = [50; -30];    % [ΔF3; ΔF4]
d_dev_B(:,k_step:end) = repmat(d_step_val,1,Tsim-k_step+1);

%% ========================================================================
% 8) Linear-plant simulations (constrained MPC)
% ========================================================================
disp('=== PROBLEM 9 – EXPERIMENT A: constrained MPC, linear plant ===');
simA = simulateOffsetFreeMPC_constrained( ...
    A,B,E,C, A_e,B_e,C_e, K_e, ...
    H,Phi,Gamma,Qbar, u_s,zs, ...
    z_ref_dev_A, d_dev_A, R_meas, lb,ub, Tsim);

disp('=== PROBLEM 9 – EXPERIMENT B: constrained MPC, linear plant ===');
simB = simulateOffsetFreeMPC_constrained( ...
    A,B,E,C, A_e,B_e,C_e, K_e, ...
    H,Phi,Gamma,Qbar, u_s,zs, ...
    z_ref_dev_B, d_dev_B, R_meas, lb,ub, Tsim);

%% ========================================================================
% 9) SDE plant overlay (same controller: constrained MPC + linear KF)
% ========================================================================
sigma3 = 2; sigma4 = 2;
sigma_meas_4 = [sqrt(R_meas(1,1)); sqrt(R_meas(2,2)); 0; 0];
x0 = xs;

disp('=== PROBLEM 9 – EXPERIMENT A: constrained MPC, SDE plant ===');
simA_sde = simulateOffsetFreeMPC_constrained_SDEplant( ...
    t, x0, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
    A_e,B_e,C_e, K_e, H,Phi,Gamma,Qbar, ...
    u_s, d_s, zs, z_ref_dev_A, d_dev_A, lb,ub, Tsim);

disp('=== PROBLEM 9 – EXPERIMENT B: constrained MPC, SDE plant ===');
simB_sde = simulateOffsetFreeMPC_constrained_SDEplant( ...
    t, x0, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
    A_e,B_e,C_e, K_e, H,Phi,Gamma,Qbar, ...
    u_s, d_s, zs, z_ref_dev_B, d_dev_B, lb,ub, Tsim);

%% ========================================================================
% 10) PLOTS (like your Problem 8 style + step marker)
% ========================================================================

% -------- Exp A outputs --------
figure('Name','P9_ExpA_outputs');
sgtitle('Problem 9 – Experiment A: reference steps, d\_dev=0 (linear vs SDE)','Interpreter','none');

subplot(3,1,1); hold on; grid on;
plot(t, z_ref_abs_A(1,:), 'k--','LineWidth',1.2,'DisplayName','h1 ref');
plot(t, simA.y_abs(1,:),  'b', 'LineWidth',1.4,'DisplayName','h1 linear');
plot(t, simA_sde.y_abs(1,:), 'm:', 'LineWidth',1.4,'DisplayName','h1 SDE meas');
ylabel('$h_1$ [cm]','Interpreter','latex');
legend('Location','best');

subplot(3,1,2); hold on; grid on;
plot(t, z_ref_abs_A(2,:), 'k--','LineWidth',1.2,'DisplayName','h2 ref');
plot(t, simA.y_abs(2,:),  'b', 'LineWidth',1.4,'DisplayName','h2 linear');
plot(t, simA_sde.y_abs(2,:), 'm:', 'LineWidth',1.4,'DisplayName','h2 SDE meas');
ylabel('$h_2$ [cm]','Interpreter','latex');
legend('Location','best');

subplot(3,1,3); hold on; grid on;
plot(t, simA.u(1,:), 'b','LineWidth',1.4,'DisplayName','u1 linear');
plot(t, simA.u(2,:), 'r','LineWidth',1.4,'DisplayName','u2 linear');
plot(t, simA_sde.u(1,:), 'm:','LineWidth',1.2,'DisplayName','u1 SDE');
plot(t, simA_sde.u(2,:), 'm-.','LineWidth',1.2,'DisplayName','u2 SDE');
yline(u_s(1),'k--','LineWidth',1.0,'HandleVisibility','off');
yline(u_s(2),'k--','LineWidth',1.0,'HandleVisibility','off');
yline(0,'k:','LineWidth',1.0,'HandleVisibility','off');
yline(600,'k:','LineWidth',1.0,'HandleVisibility','off');
ylabel('$u$ [cm$^3$/s]','Interpreter','latex');
xlabel('Time [s]');
legend('Location','best');

% -------- Exp A disturbance estimates --------
figure('Name','P9_ExpA_disturbances');
sgtitle('Problem 9 – Experiment A: disturbance estimates (linear vs SDE)','Interpreter','none');

subplot(2,1,1); hold on; grid on;
plot(t, simA.d_hat(1,:), 'b','LineWidth',1.4,'DisplayName','d1 hat linear');
plot(t, simA_sde.d_hat(1,:), 'm:','LineWidth',1.2,'DisplayName','d1 hat SDE');
ylabel('$\hat d_1$ [cm$^3$/s]','Interpreter','latex');
ylim([-50 50])
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(t, simA.d_hat(2,:), 'b','LineWidth',1.4,'DisplayName','d2 hat linear');
plot(t, simA_sde.d_hat(2,:), 'm:','LineWidth',1.2,'DisplayName','d2 hat SDE');
ylabel('$\hat d_2$ [cm$^3$/s]','Interpreter','latex');
xlabel('Time [s]');
ylim([-50 50])
legend('Location','best');

% -------- Exp B outputs --------
figure('Name','P9_ExpB_outputs');
sgtitle('Problem 9 – Experiment B: stepped refs + disturbance step (linear vs SDE)','Interpreter','none');

subplot(3,1,1); hold on; grid on;
plot(t, z_ref_abs_B(1,:), 'k--','LineWidth',1.2,'DisplayName','h1 ref');
plot(t, simB.y_abs(1,:),  'b', 'LineWidth',1.4,'DisplayName','h1 linear');
plot(t, simB_sde.y_abs(1,:), 'm:', 'LineWidth',1.4,'DisplayName','h1 SDE meas');
addStepMarker(gca, t(k_step), 'step in d', 0.01, 0.60);
ylabel('$h_1$ [cm]','Interpreter','latex');
legend('Location','best');

subplot(3,1,2); hold on; grid on;
plot(t, z_ref_abs_B(2,:), 'k--','LineWidth',1.2,'DisplayName','h2 ref');
plot(t, simB.y_abs(2,:),  'b', 'LineWidth',1.4,'DisplayName','h2 linear');
plot(t, simB_sde.y_abs(2,:), 'm:', 'LineWidth',1.4,'DisplayName','h2 SDE meas');
addStepMarker(gca, t(k_step), 'step in d', 0.01, 0.60);
ylabel('$h_2$ [cm]','Interpreter','latex');
legend('Location','best');

subplot(3,1,3); hold on; grid on;
plot(t, simB.u(1,:), 'b','LineWidth',1.4,'DisplayName','u1 linear');
plot(t, simB.u(2,:), 'r','LineWidth',1.4,'DisplayName','u2 linear');
plot(t, simB_sde.u(1,:), 'm:','LineWidth',1.2,'DisplayName','u1 SDE');
plot(t, simB_sde.u(2,:), 'm-.','LineWidth',1.2,'DisplayName','u2 SDE');
addStepMarker(gca, t(k_step), 'step in d', 0.01, 0.60);
yline(u_s(1),'k--','LineWidth',1.0,'HandleVisibility','off');
yline(u_s(2),'k--','LineWidth',1.0,'HandleVisibility','off');
yline(0,'k:','LineWidth',1.0,'HandleVisibility','off');
yline(600,'k:','LineWidth',1.0,'HandleVisibility','off');
ylabel('$u$ [cm$^3$/s]','Interpreter','latex');
xlabel('Time [s]');
legend('Location','best');

% -------- Exp B disturbances --------
figure('Name','P9_ExpB_disturbances');
sgtitle('Problem 9 – Experiment B: true vs estimated disturbances (linear vs SDE)','Interpreter','none');

subplot(2,1,1); hold on; grid on;
plot(t, d_dev_B(1,:), 'k--','LineWidth',1.4,'DisplayName','d1 true');
plot(t, simB.d_hat(1,:), 'b','LineWidth',1.4,'DisplayName','d1 hat linear');
plot(t, simB_sde.d_hat(1,:), 'm:','LineWidth',1.2,'DisplayName','d1 hat SDE');
addStepMarker(gca, t(k_step), 'step in d', 0.01, 0.60);
ylabel('$\Delta F_3$ [cm$^3$/s]','Interpreter','latex');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(t, d_dev_B(2,:), 'k--','LineWidth',1.4,'DisplayName','d2 true');
plot(t, simB.d_hat(2,:), 'b','LineWidth',1.4,'DisplayName','d2 hat linear');
plot(t, simB_sde.d_hat(2,:), 'm:','LineWidth',1.2,'DisplayName','d2 hat SDE');
addStepMarker(gca, t(k_step), 'step in d', 0.01, 0.60);
ylabel('$\Delta F_4$ [cm$^3$/s]','Interpreter','latex');
xlabel('Time [s]');
legend('Location','best');

%% ========================================================================
% 11) Save all figures as PDFs
% ========================================================================
saveAllFiguresAsPDF(figDir);
fprintf('Saved figures to: %s\n', figDir);

%% ========================================================================
% ===================== HELPER FUNCTIONS ==================================
% ========================================================================

function [Phi_x, Gamma_x] = buildPredictionMatrices(A,B,N)
    [nx,nu] = size(B);
    Phi_x   = zeros(N*nx, nx);
    Gamma_x = zeros(N*nx, N*nu);

    A_power = eye(nx);
    for i = 1:N
        A_power = A * A_power;          % A^i
        rows = (i-1)*nx + (1:nx);
        Phi_x(rows,:) = A_power;

        A_j = eye(nx);
        for j = 1:i
            cols = (j-1)*nu + (1:nu);
            Gamma_x(rows,cols) = A_j * B;
            A_j = A * A_j;
        end
    end
end

function [h1_ref, h2_ref] = generate_stair_references(t, t_change_h1, values_h1, t_change_h2, values_h2)
    Nt = numel(t);
    h1_ref = zeros(1,Nt);
    h2_ref = zeros(1,Nt);

    for k = 1:Nt
        tk = t(k);

        idx1 = find(tk >= t_change_h1, 1, 'last');
        if isempty(idx1), idx1 = 1; end
        idx1 = min(idx1, numel(values_h1));
        h1_ref(k) = values_h1(idx1);

        idx2 = find(tk >= t_change_h2, 1, 'last');
        if isempty(idx2), idx2 = 1; end
        idx2 = min(idx2, numel(values_h2));
        h2_ref(k) = values_h2(idx2);
    end
end

function addStepMarker(ax, t_step, labelStr, dxFrac, yFromTopFrac)
% Draw vertical red dashed line at t_step, and put text near it.
% dxFrac: fraction of x-range to shift text to the right (use negative to go left)
% yFromTopFrac: fraction of y-range DOWN from top (0.1 near top, 0.6 mid)
    if nargin < 4, dxFrac = 0.01; end
    if nargin < 5, yFromTopFrac = 0.60; end

    xl = xlim(ax);
    yl = ylim(ax);

    line(ax, [t_step t_step], yl, ...
        'LineStyle','--','LineWidth',1.2, 'Color',[1 0 0], ...
        'HandleVisibility','off');

    x_text = t_step + dxFrac*(xl(2)-xl(1));
    y_text = yl(2) - yFromTopFrac*(yl(2)-yl(1));

    text(ax, x_text, y_text, labelStr, ...
        'Color',[1 0 0], 'FontSize',12, ...
        'Interpreter','none', 'Clipping','on');
end

function sim = simulateOffsetFreeMPC_constrained( ...
            A,B,E,C, A_e,B_e,C_e, K_e, ...
            H,Phi,Gamma,Qbar, u_s,zs, ...
            z_ref_dev, d_dev, R_meas, lb,ub, Tsim)
% Constrained offset-free MPC + static KF on *linear discrete plant*.
    [nx,nu] = size(B);
    nd      = size(E,2);
    ny      = size(C,1);
    nx_e    = nx + nd;
    N       = size(Phi,1)/ny;

    x_dev   = zeros(nx, Tsim+1);
    y_true  = zeros(ny, Tsim);
    y_abs   = zeros(ny, Tsim);
    y_meas  = zeros(ny, Tsim);

    x_hat_e = zeros(nx_e, Tsim+1);
    u_dev   = zeros(nu, Tsim);
    u_abs   = zeros(nu, Tsim);

    sigma_y = sqrt(diag(R_meas));

    % QP parts (bounds only)
    l = lb; u = ub;
    Aqp = []; bl = []; bu = [];
    xinit = zeros(N*nu,1);

    for k = 1:Tsim
        if k > 1
            x_dev(:,k) = A*x_dev(:,k-1) + B*u_dev(:,k-1) + E*d_dev(:,k-1);
        end

        y_true(:,k) = C*x_dev(:,k);
        y_abs(:,k)  = zs + y_true(:,k);
        y_meas(:,k) = y_true(:,k) + sigma_y .* randn(ny,1);

        % KF meas update
        innov        = y_meas(:,k) - C_e*x_hat_e(:,k);
        x_hat_e(:,k) = x_hat_e(:,k) + K_e*innov;

        % stack reference
        zref_stack = zeros(ny*N,1);
        for j = 0:N-1
            idx = min(k+j, Tsim);
            zref_stack(j*ny+(1:ny)) = z_ref_dev(:,idx);
        end

        % f for 0.5 U'HU + f'U
        f = 2*Gamma'*Qbar*(Phi*x_hat_e(:,k) - zref_stack);

        [U_seq, ~] = problem_7(H, f, l, u, Aqp, bl, bu, xinit);
        u_dev(:,k) = U_seq(1:nu);
        u_abs(:,k) = u_s + u_dev(:,k);

        % KF time update
        if k < Tsim
            x_hat_e(:,k+1) = A_e*x_hat_e(:,k) + B_e*u_dev(:,k);
        end
    end

    sim = struct();
    sim.y_abs = y_abs;
    sim.u     = u_abs;
    sim.d_hat = x_hat_e(nx+1:end,1:Tsim);
end

function sim = simulateOffsetFreeMPC_constrained_SDEplant( ...
            t, x0, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
            A_e,B_e,C_e, K_e, H,Phi,Gamma,Qbar, ...
            u_s, d_s, zs, z_ref_dev, d_dev, lb,ub, Tsim)
% Constrained offset-free MPC + static KF (linear) on NONLINEAR SDE plant.
    Ts = t(2) - t(1);
    ny = length(output_index);
    nu = 2;
    nd = 2;
    nx_e = 4 + nd;
    N = size(Phi,1)/ny;

    tf = t(end);
    [~,~,dW3] = ScalarStdWienerProcess(tf, length(t), 1);
    [~,~,dW4] = ScalarStdWienerProcess(tf, length(t), 1);
    if length(dW3) < Tsim || length(dW4) < Tsim
        error('ScalarStdWienerProcess returned too few increments.');
    end

    x_true = zeros(4, Tsim+1);
    y_meas_abs = zeros(ny, Tsim);

    x_hat_e = zeros(nx_e, Tsim+1);
    u_abs = zeros(nu, Tsim);
    u_dev = zeros(nu, Tsim);

    x_true(:,1) = x0;

    % QP parts (bounds only)
    l = lb; u = ub;
    Aqp = []; bl = []; bu = [];
    xinit = zeros(N*nu,1);

    for k = 1:Tsim
        % measure plant
        y_true_abs4 = FourTankSystemSensor(x_true(:,k), p);
        y_meas_abs(:,k) = y_true_abs4(output_index) + sigma_meas_4(output_index).*randn(ny,1);
        y_meas_dev = y_meas_abs(:,k) - zs;

        % KF meas update
        innov        = y_meas_dev - C_e*x_hat_e(:,k);
        x_hat_e(:,k) = x_hat_e(:,k) + K_e*innov;

        % stack reference
        zref_stack = zeros(ny*N,1);
        for j = 0:N-1
            idx = min(k+j, Tsim);
            zref_stack(j*ny+(1:ny)) = z_ref_dev(:,idx);
        end

        % QP
        f = 2*Gamma'*Qbar*(Phi*x_hat_e(:,k) - zref_stack);
        [U_seq, ~] = problem_7(H, f, l, u, Aqp, bl, bu, xinit);

        u_dev(:,k) = U_seq(1:nu);
        u_abs(:,k) = u_s + u_dev(:,k);

        % KF time update
        if k < Tsim
            x_hat_e(:,k+1) = A_e*x_hat_e(:,k) + B_e*u_dev(:,k);
        end

        % plant step (Euler–Maruyama)
        if k < Tsim
            d_abs_k = d_s + d_dev(:,k);
            f_nl = FourTankSystemModified(0, x_true(:,k), u_abs(:,k), d_abs_k, p);
            noise = [0; 0; sigma3*dW3(k); sigma4*dW4(k)];
            x_true(:,k+1) = x_true(:,k) + Ts*f_nl + noise;
        end
    end

    sim = struct();
    sim.y_abs = y_meas_abs;                 % measured ABS heights (subset)
    sim.u     = u_abs;
    sim.d_hat = x_hat_e(4+1:end,1:Tsim);
end

function saveAllFiguresAsPDF(figDir)
    figs = findall(0,'Type','figure');
    for i = 1:numel(figs)
        fig = figs(i);
        figName = get(fig,'Name');
        if isempty(figName)
            figName = sprintf('Figure_%d', fig.Number);
        end
        figName = regexprep(figName,'[^a-zA-Z0-9_\-]','_');
        outFile = fullfile(figDir, [figName '.pdf']);

        try
            exportgraphics(fig, outFile, 'ContentType','vector');
        catch
            print(fig, outFile, '-dpdf', '-painters');
        end
    end
end
