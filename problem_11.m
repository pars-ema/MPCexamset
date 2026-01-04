%% ========================================================================
% compare_P8_P9_P10.m
% Fair comparison of offset-free MPCs:
%   P8  = unconstrained MPC (analytic solution)
%   P9  = hard input bounds (QP with bounds)
%   P10 = hard input bounds + hard move bounds + SOFT output upper bound (slacks)
%
% Runs:
%   Experiment A: reference steps, d_dev = 0
%   Experiment B: same reference steps + step in disturbance (F3,F4)
%
% NEW (requested):
%   - Stores BOTH true heights and measured heights for ALL sims:
%       sim.y_true_abs (ny x Tsim)
%       sim.y_meas_abs (ny x Tsim)
%   - Plots true as lines + measured as sparse 'x' markers.
%
% Produces and saves vector-PDF figures in:
%   figures/compare_8_9_10/
%
% REQUIREMENTS ON PATH:
%   - FourTankSystemModified.m
%   - FourTankSystemSensor.m
%   - ScalarStdWienerProcess.m
%   - problem_7.m     (your QP solver interface wrapping quadprog)
% ========================================================================

clear; clc; close all;

%% ------------------------- output folder + style -------------------------
figDir = fullfile('figures','compare_8_9_10');
if ~exist(figDir,'dir'), mkdir(figDir); end

set(groot,'defaultFigureColor','w');
set(groot,'defaultFigurePaperPositionMode','auto');
set(groot,'defaultFigurePosition',[100 100 1400 850]);

set(groot,'defaultTextInterpreter','none');
set(groot,'defaultLegendInterpreter','none');
set(groot,'defaultAxesTickLabelInterpreter','none');
set(groot,'defaultAxesFontSize',13);
set(groot,'defaultLegendFontSize',12);

rng(0);  % reproducibility for measurement noise

%% ========================================================================
% 0) Common physical parameters + steady state
% ========================================================================
Ts = 10;                           % sampling [s]
nx = 4; nu = 2; nd = 2; ny = 2;

a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327;
g  = 981;
rho = 1;
gamma1 = 0.58;
gamma2 = 0.72;

p = [a1;a2;a3;a4;A1;A2;A3;A4;g;gamma1;gamma2;rho];

u_s = [300;300];                   % steady input (F1,F2)
d_s = [250;250];                   % steady disturbance (F3,F4) for NL plant
output_index = [1 2];

xs_guess = 500*ones(4,1);
opts = optimoptions('fsolve','Display','none');
xs = fsolve(@(x)FourTankSystemModified(0,x,u_s,d_s,p), xs_guess, opts);

h_s_abs_all = FourTankSystemSensor(xs,p);
zs = h_s_abs_all(output_index);    % measured steady outputs (h1,h2)

fprintf("Steady-state xs:\n"); disp(xs);
fprintf("Steady-state zs=[h1;h2]:\n"); disp(zs);

%% ========================================================================
% 1) Linearization + ZOH discretization (common to all)
% ========================================================================
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

A  = Md(1:nx,1:nx);
BEd = Md(1:nx,nx+1:end);
B  = BEd(:,1:nu);
E  = BEd(:,nu+1:end);
C  = C_c;

%% ========================================================================
% 2) Augmentation (offset-free) + static KF (common)
% ========================================================================
Qw     = diag([25 25]);
Qx     = E * Qw * E';
Qd_RW  = 1 * eye(nd);
R_meas = diag([4 4]);

A_e = [A, E;
       zeros(nd,nx), eye(nd)];
B_e = [B;
       zeros(nd,nu)];
C_e = [C, zeros(ny,nd)];

Q_e = blkdiag(Qx, Qd_RW);
K_e = computeKF(A_e,C_e,Q_e,R_meas);
fprintf("Static augmented KF gain norm: %.4f\n", norm(K_e));

%% ========================================================================
% 3) Common MPC settings (same horizon/weights)
% ========================================================================
N = 15;

Wz  = diag([1 1]);
Wu  = 1e-4*eye(nu);
Wdu = 1e-7*eye(nu);          % only used by P10 (move penalty)

% prediction matrices for augmented model
[Phi_x_e, Phi_u_e] = predMatrices(A_e,B_e,N);
Phi   = kron(eye(N),C_e) * Phi_x_e;
Gamma = kron(eye(N),C_e) * Phi_u_e;

Qbar = kron(eye(N),Wz);
Rbar = kron(eye(N),Wu);

% P8/P9 quadratic term (outputs+inputs)
H_u = 2*(Gamma'*Qbar*Gamma + Rbar);

%% ========================================================================
% 4) Experiments A/B definition (COMMON)
% ========================================================================
Tsim = 300;
t    = (0:Tsim-1)*Ts;

% Reference schedule (use ONE schedule for all controllers)
t_change_h1 = [0 600 1600];
values_h1   = [zs(1) zs(1)*1.3 zs(1)*0.8];

t_change_h2 = [0 600 1300];
values_h2   = [zs(2) zs(2)*1.3 zs(2)*1.0];

[h1_ref, h2_ref] = generate_stair_references(t, t_change_h1, values_h1, ...
                                                t_change_h2, values_h2);
z_ref_abs = [h1_ref; h2_ref];
z_ref_dev = z_ref_abs - zs;

% Disturbances (deviation form)
d_dev_A = zeros(nd,Tsim);

d_dev_B = zeros(nd,Tsim);
k_step  = round(2000/Ts);
k_step  = max(1, min(Tsim, k_step));
d_step_val = [50; -20];                  % [ΔF3; ΔF4]
d_dev_B(:,k_step:end) = repmat(d_step_val,1,Tsim-k_step+1);

% SDE noise settings (nonlinear plant)
sigma3 = 2; sigma4 = 2;
sigma_meas_4 = [sqrt(R_meas(1,1)); sqrt(R_meas(2,2)); 0; 0];

%% ========================================================================
% 5) Constraints (use ONE consistent set for P9 + P10)
% ========================================================================
% Hard input amplitude bounds (absolute pumps)
u_min_abs = [0;0];
u_max_abs = [600;600];

lb_u_dev = u_min_abs - u_s;
ub_u_dev = u_max_abs - u_s;

lb_stack = repmat(lb_u_dev, N, 1);
ub_stack = repmat(ub_u_dev, N, 1);

% Hard move bounds for P10 only (absolute move in deviation space)
du_max_abs = [300;300];
du_min_abs = -du_max_abs;

du_min_stack = kron(ones(N,1),du_min_abs);
du_max_stack = kron(ones(N,1),du_max_abs);

% Soft output upper bounds (absolute) for P10
z_max_abs = [150;150];
z_max_dev = z_max_abs - zs;
Z_max_stack = kron(ones(N,1), z_max_dev);

%% ========================================================================
% 6) RUN ALL CONTROLLERS — Experiment A (linear + SDE)
% ========================================================================
disp('=== Experiment A: Linear ===');
simA_P8_lin  = simP8_linear(A,B,E,C, A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H_u, u_s,zs, z_ref_dev, d_dev_A, R_meas, Tsim);
simA_P9_lin  = simP9_linear(A,B,E,C, A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H_u, u_s,zs, z_ref_dev, d_dev_A, R_meas, lb_stack,ub_stack, Tsim);
simA_P10_lin = simP10_linear(A,B,E,C, A_e,B_e,C_e, K_e, Phi_x_e,Phi_u_e, Wz,Wu,Wdu, ...
                             u_s,zs, z_ref_dev, d_dev_A, R_meas, ...
                             lb_stack,ub_stack, du_min_stack,du_max_stack, Z_max_stack, Tsim);

disp('=== Experiment A: SDE plant ===');
simA_P8_sde  = simP8_SDE(t, xs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
                        A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H_u, u_s,d_s,zs, z_ref_dev, d_dev_A, Tsim);
simA_P9_sde  = simP9_SDE(t, xs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
                        A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H_u, u_s,d_s,zs, z_ref_dev, d_dev_A, lb_stack,ub_stack, Tsim);
simA_P10_sde = simP10_SDE(t, xs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
                          A,B,E,C, A_e,B_e,C_e, K_e, Phi_x_e,Phi_u_e, Wz,Wu,Wdu, ...
                          u_s,d_s,zs, z_ref_dev, d_dev_A, ...
                          lb_stack,ub_stack, du_min_stack,du_max_stack, Z_max_stack, Tsim);

%% ========================================================================
% 7) RUN ALL CONTROLLERS — Experiment B (linear + SDE)
% ========================================================================
disp('=== Experiment B: Linear ===');
simB_P8_lin  = simP8_linear(A,B,E,C, A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H_u, u_s,zs, z_ref_dev, d_dev_B, R_meas, Tsim);
simB_P9_lin  = simP9_linear(A,B,E,C, A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H_u, u_s,zs, z_ref_dev, d_dev_B, R_meas, lb_stack,ub_stack, Tsim);
simB_P10_lin = simP10_linear(A,B,E,C, A_e,B_e,C_e, K_e, Phi_x_e,Phi_u_e, Wz,Wu,Wdu, ...
                             u_s,zs, z_ref_dev, d_dev_B, R_meas, ...
                             lb_stack,ub_stack, du_min_stack,du_max_stack, Z_max_stack, Tsim);

disp('=== Experiment B: SDE plant ===');
simB_P8_sde  = simP8_SDE(t, xs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
                        A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H_u, u_s,d_s,zs, z_ref_dev, d_dev_B, Tsim);
simB_P9_sde  = simP9_SDE(t, xs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
                        A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H_u, u_s,d_s,zs, z_ref_dev, d_dev_B, lb_stack,ub_stack, Tsim);
simB_P10_sde = simP10_SDE(t, xs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
                          A,B,E,C, A_e,B_e,C_e, K_e, Phi_x_e,Phi_u_e, Wz,Wu,Wdu, ...
                          u_s,d_s,zs, z_ref_dev, d_dev_B, ...
                          lb_stack,ub_stack, du_min_stack,du_max_stack, Z_max_stack, Tsim);

%% ========================================================================
% 8) PLOTS (comparison overlays) + save PDFs
% ========================================================================

% Marker thinning for measured crosses (to avoid clutter)
mkStep = 4;   % plot every mkStep samples as 'x'

% -------- Exp A: linear comparison --------
plotCompareOutputsInputs(t, z_ref_abs, simA_P8_lin, simA_P9_lin, simA_P10_lin, ...
    u_s, u_min_abs, u_max_abs, [], ...
    'ExpA_LINEAR__P8_vs_P9_vs_P10', ...
    'Experiment A (Linear): P8 vs P9 vs P10', mkStep);

% -------- Exp A: SDE comparison --------
plotCompareOutputsInputs(t, z_ref_abs, simA_P8_sde, simA_P9_sde, simA_P10_sde, ...
    u_s, u_min_abs, u_max_abs, [], ...
    'ExpA_SDE__P8_vs_P9_vs_P10', ...
    'Experiment A (SDE plant): P8 vs P9 vs P10', mkStep);

% -------- Exp B: linear comparison (with step marker) --------
plotCompareOutputsInputs(t, z_ref_abs, simB_P8_lin, simB_P9_lin, simB_P10_lin, ...
    u_s, u_min_abs, u_max_abs, t(k_step), ...
    'ExpB_LINEAR__P8_vs_P9_vs_P10', ...
    'Experiment B (Linear): P8 vs P9 vs P10', mkStep);

% -------- Exp B: SDE comparison (with step marker) --------
plotCompareOutputsInputs(t, z_ref_abs, simB_P8_sde, simB_P9_sde, simB_P10_sde, ...
    u_s, u_min_abs, u_max_abs, t(k_step), ...
    'ExpB_SDE__P8_vs_P9_vs_P10', ...
    'Experiment B (SDE plant): P8 vs P9 vs P10', mkStep);

% -------- Exp B: disturbance estimates (linear) --------
plotCompareDisturbances(t, d_dev_B, simB_P8_lin, simB_P9_lin, simB_P10_lin, ...
    t(k_step), 'ExpB_LINEAR_disturbance_estimates', ...
    'Experiment B (Linear): true vs estimated disturbances');

% -------- Exp B: disturbance estimates (SDE) --------
plotCompareDisturbances(t, d_dev_B, simB_P8_sde, simB_P9_sde, simB_P10_sde, ...
    t(k_step), 'ExpB_SDE_disturbance_estimates', ...
    'Experiment B (SDE plant): true vs estimated disturbances');

% -------- P10 slack activity (Exp B recommended) --------
plotP10Slack(t, simB_P10_lin, simB_P10_sde, t(k_step), 'ExpB_P10_slack_activity', ...
    'P10 slack activity (Experiment B)');

% Save all figures
saveAllFiguresAsPDF(figDir);
fprintf('Saved comparison figures to: %s\n', figDir);

%% ========================================================================
% ============================ HELPERS ====================================
% ========================================================================

function K = computeKF(A,C,Q,R)
    [P,~,~] = dare(A',C',Q,R);
    P = 0.5*(P+P');
    Re = C*P*C' + R;
    K  = P*C'/Re;
end

function [Phi_x, Phi_u] = predMatrices(A,B,N)
    % x_{k+i} = Phi_x(i)*x_k + sum_j Phi_u(i,j)*u_{k+j}
    [nx,nu] = size(B);
    Phi_x = zeros(N*nx, nx);
    Phi_u = zeros(N*nx, N*nu);

    A_power = eye(nx);
    for i = 1:N
        A_power = A*A_power;
        rr = (i-1)*nx + (1:nx);
        Phi_x(rr,:) = A_power;

        A_j = eye(nx);
        for j = 1:i
            cc = (j-1)*nu + (1:nu);
            Phi_u(rr,cc) = A_j*B;
            A_j = A*A_j;
        end
    end
end

function [h1_ref, h2_ref] = generate_stair_references(t, t_change_h1, values_h1, t_change_h2, values_h2)
    Nt = numel(t);
    h1_ref = zeros(1,Nt);
    h2_ref = zeros(1,Nt);
    for k=1:Nt
        idx1 = find(t(k) >= t_change_h1, 1, 'last'); if isempty(idx1), idx1 = 1; end
        idx2 = find(t(k) >= t_change_h2, 1, 'last'); if isempty(idx2), idx2 = 1; end
        idx1 = min(idx1,numel(values_h1));
        idx2 = min(idx2,numel(values_h2));
        h1_ref(k) = values_h1(idx1);
        h2_ref(k) = values_h2(idx2);
    end
end

%% ------------------------- P8 (unconstrained) ---------------------------

function sim = simP8_linear(A,B,E,C, A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H, u_s,zs, z_ref_dev, d_dev, R_meas, Tsim)
    [nx,nu] = size(B);
    nd = size(E,2);
    ny = size(C,1);
    nx_e = nx+nd;
    N = size(Phi,1)/ny;

    x_dev   = zeros(nx, Tsim+1);
    y_abs   = zeros(ny, Tsim);
    u_abs   = zeros(nu, Tsim);
    x_hat_e = zeros(nx_e, Tsim+1);

    % NEW:
    y_true_abs = zeros(ny,Tsim);
    y_meas_abs = zeros(ny,Tsim);

    sigma_y = sqrt(diag(R_meas));

    for k=1:Tsim
        if k>1
            x_dev(:,k) = A*x_dev(:,k-1) + B*(u_abs(:,k-1)-u_s) + E*d_dev(:,k-1);
        end
        y_dev = C*x_dev(:,k);
        y_abs(:,k) = zs + y_dev;

        y_meas = y_dev + sigma_y.*randn(ny,1);   % deviation measurement

        % NEW store true/meas in ABS units
        y_true_abs(:,k) = y_abs(:,k);
        y_meas_abs(:,k) = zs + y_meas;

        % KF update
        innov = y_meas - C_e*x_hat_e(:,k);
        x_hat_e(:,k) = x_hat_e(:,k) + K_e*innov;

        % stack reference
        zref_stack = stackRef(z_ref_dev, k, N);
        g = 2*Gamma'*Qbar*(Phi*x_hat_e(:,k) - zref_stack);

        U = -H \ g;                  % unconstrained
        u_dev_k = U(1:nu);
        u_abs(:,k) = u_s + u_dev_k;

        if k<Tsim
            x_hat_e(:,k+1) = A_e*x_hat_e(:,k) + B_e*u_dev_k;
        end
    end

    sim = struct();
    sim.y_abs      = y_abs;           % kept (for legacy)
    sim.y_true_abs = y_true_abs;      % NEW
    sim.y_meas_abs = y_meas_abs;      % NEW
    sim.u_abs      = u_abs;
    sim.d_hat      = x_hat_e(nx+1:end,1:Tsim);
    sim.eps        = zeros(ny,Tsim);
end

function sim = simP8_SDE(t, x0_abs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
                         A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H, u_s,d_s,zs, z_ref_dev, d_dev, Tsim)
    ny = numel(output_index);
    nu = 2; nd = 2;
    nx_e = size(A_e,1);
    N = size(Phi,1)/ny;

    tf = t(end);
    [~,~,dW3] = ScalarStdWienerProcess(tf, numel(t), 1);
    [~,~,dW4] = ScalarStdWienerProcess(tf, numel(t), 1);

    x_true = zeros(4,Tsim+1);
    x_true(:,1) = x0_abs;

    y_abs = zeros(ny,Tsim);
    u_abs = zeros(nu,Tsim);
    x_hat_e = zeros(nx_e,Tsim+1);

    % NEW:
    y_true_abs = zeros(ny,Tsim);
    y_meas_abs_store = zeros(ny,Tsim);

    for k=1:Tsim
        y_true_abs4 = FourTankSystemSensor(x_true(:,k), p);
        y_meas_abs  = y_true_abs4(output_index) + sigma_meas_4(output_index).*randn(ny,1);
        y_abs(:,k)  = y_meas_abs;   % measured shown in plot usually

        % NEW store:
        y_true_abs(:,k) = y_true_abs4(output_index);
        y_meas_abs_store(:,k) = y_meas_abs;

        y_meas_dev  = y_meas_abs - zs;

        innov = y_meas_dev - C_e*x_hat_e(:,k);
        x_hat_e(:,k) = x_hat_e(:,k) + K_e*innov;

        zref_stack = stackRef(z_ref_dev, k, N);
        g = 2*Gamma'*Qbar*(Phi*x_hat_e(:,k) - zref_stack);
        U = -H \ g;

        u_dev_k = U(1:nu);
        u_k_abs = u_s + u_dev_k;
        u_abs(:,k) = u_k_abs;

        if k<Tsim
            x_hat_e(:,k+1) = A_e*x_hat_e(:,k) + B_e*u_dev_k;

            d_abs_k = d_s + d_dev(:,k);
            f = FourTankSystemModified(0, x_true(:,k), u_k_abs, d_abs_k, p);
            noise = [0;0; sigma3*dW3(k); sigma4*dW4(k)];
            x_true(:,k+1) = x_true(:,k) + (t(2)-t(1))*f + noise;
        end
    end

    sim = struct();
    sim.y_abs      = y_abs;
    sim.y_true_abs = y_true_abs;          % NEW
    sim.y_meas_abs = y_meas_abs_store;    % NEW
    sim.u_abs      = u_abs;
    sim.d_hat      = x_hat_e(end-nd+1:end,1:Tsim);
    sim.eps        = zeros(ny,Tsim);
end

%% ------------------------- P9 (hard input bounds) -----------------------

function sim = simP9_linear(A,B,E,C, A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H, u_s,zs, z_ref_dev, d_dev, R_meas, lb,ub, Tsim)
    [nx,nu] = size(B);
    nd = size(E,2);
    ny = size(C,1);
    nx_e = nx+nd;
    N = size(Phi,1)/ny;

    x_dev   = zeros(nx, Tsim+1);
    y_abs   = zeros(ny, Tsim);
    u_abs   = zeros(nu, Tsim);
    x_hat_e = zeros(nx_e, Tsim+1);

    % NEW:
    y_true_abs = zeros(ny,Tsim);
    y_meas_abs = zeros(ny,Tsim);

    sigma_y = sqrt(diag(R_meas));

    Aqp = []; bl = []; bu = [];
    xinit = zeros(N*nu,1);

    for k=1:Tsim
        if k>1
            x_dev(:,k) = A*x_dev(:,k-1) + B*(u_abs(:,k-1)-u_s) + E*d_dev(:,k-1);
        end
        y_dev = C*x_dev(:,k);
        y_abs(:,k) = zs + y_dev;

        y_meas = y_dev + sigma_y.*randn(ny,1);

        % NEW store true/meas ABS
        y_true_abs(:,k) = y_abs(:,k);
        y_meas_abs(:,k) = zs + y_meas;

        innov = y_meas - C_e*x_hat_e(:,k);
        x_hat_e(:,k) = x_hat_e(:,k) + K_e*innov;

        zref_stack = stackRef(z_ref_dev, k, N);
        f = 2*Gamma'*Qbar*(Phi*x_hat_e(:,k) - zref_stack);

        [U,~] = problem_7(H, f, lb, ub, Aqp, bl, bu, xinit);
        u_dev_k = U(1:nu);
        u_abs(:,k) = u_s + u_dev_k;

        if k<Tsim
            x_hat_e(:,k+1) = A_e*x_hat_e(:,k) + B_e*u_dev_k;
        end
    end

    sim = struct();
    sim.y_abs      = y_abs;
    sim.y_true_abs = y_true_abs;   % NEW
    sim.y_meas_abs = y_meas_abs;   % NEW
    sim.u_abs      = u_abs;
    sim.d_hat      = x_hat_e(end-nd+1:end,1:Tsim);
    sim.eps        = zeros(ny,Tsim);
end

function sim = simP9_SDE(t, x0_abs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
                         A_e,B_e,C_e, K_e, Phi,Gamma,Qbar,H, u_s,d_s,zs, z_ref_dev, d_dev, lb,ub, Tsim)
    %#ok<INUSD>
    ny = numel(output_index);
    nu = 2; nd = 2;
    nx_e = size(A_e,1);
    N = size(Phi,1)/ny;

    tf = t(end);
    [~,~,dW3] = ScalarStdWienerProcess(tf, numel(t), 1);
    [~,~,dW4] = ScalarStdWienerProcess(tf, numel(t), 1);

    x_true = zeros(4,Tsim+1);
    x_true(:,1) = x0_abs;

    y_abs = zeros(ny,Tsim);
    u_abs = zeros(nu,Tsim);
    x_hat_e = zeros(nx_e,Tsim+1);

    % NEW:
    y_true_abs = zeros(ny,Tsim);
    y_meas_abs_store = zeros(ny,Tsim);

    Aqp = []; bl = []; bu = [];
    xinit = zeros(N*nu,1);

    for k=1:Tsim
        y_true_abs4 = FourTankSystemSensor(x_true(:,k), p);
        y_meas_abs  = y_true_abs4(output_index) + sigma_meas_4(output_index).*randn(ny,1);
        y_abs(:,k)  = y_meas_abs;

        % NEW store:
        y_true_abs(:,k) = y_true_abs4(output_index);
        y_meas_abs_store(:,k) = y_meas_abs;

        y_meas_dev  = y_meas_abs - zs;

        innov = y_meas_dev - C_e*x_hat_e(:,k);
        x_hat_e(:,k) = x_hat_e(:,k) + K_e*innov;

        zref_stack = stackRef(z_ref_dev, k, N);
        f = 2*Gamma'*Qbar*(Phi*x_hat_e(:,k) - zref_stack);

        [U,~] = problem_7(H, f, lb, ub, Aqp, bl, bu, xinit);

        u_dev_k = U(1:nu);
        u_k_abs = u_s + u_dev_k;
        u_abs(:,k) = u_k_abs;

        if k<Tsim
            x_hat_e(:,k+1) = A_e*x_hat_e(:,k) + B_e*u_dev_k;

            d_abs_k = d_s + d_dev(:,k);
            fplant = FourTankSystemModified(0, x_true(:,k), u_k_abs, d_abs_k, p);
            noise = [0;0; sigma3*dW3(k); sigma4*dW4(k)];
            x_true(:,k+1) = x_true(:,k) + (t(2)-t(1))*fplant + noise;
        end
    end

    sim = struct();
    sim.y_abs      = y_abs;
    sim.y_true_abs = y_true_abs;          % NEW
    sim.y_meas_abs = y_meas_abs_store;    % NEW
    sim.u_abs      = u_abs;
    sim.d_hat      = x_hat_e(end-nd+1:end,1:Tsim);
    sim.eps        = zeros(ny,Tsim);
end

%% -------------------- P10 (u bounds + du bounds + soft zmax) ------------

function sim = simP10_linear(A,B,E,C, A_e,B_e,C_e, K_e, Phi_x_e,Phi_u_e, Wz,Wu,Wdu, ...
                             u_s,zs, z_ref_dev, d_dev, R_meas, ...
                             lb_u_stack,ub_u_stack, du_min_stack,du_max_stack, Z_max_stack, Tsim)
    [nx,nu] = size(B);
    nd = size(E,2);
    ny = size(C,1);
    nx_e = nx+nd;
    N = numel(lb_u_stack)/nu;

    % Output prediction matrices
    Phi_z_x = kron(eye(N),C_e) * Phi_x_e;
    Phi_z_u = kron(eye(N),C_e) * Phi_u_e;

    Wz_big  = kron(eye(N),Wz);
    Wu_big  = kron(eye(N),Wu);
    Wdu_big = kron(eye(N),Wdu);

    % Delta-u mapping
    E_du = makeEdU(N,nu);

    nU = N*nu;
    nZ = N*ny;
    nEta = nZ;

    % slack penalty
    w_eta = 1e3;
    Weta = w_eta*eye(nEta);

    % Hessian
    H_z  = Phi_z_u' * Wz_big  * Phi_z_u;
    H_u  = Wu_big;
    H_du = E_du'   * Wdu_big * E_du;
    H_UU = H_z + H_u + H_du;

    H = blkdiag(H_UU, Weta);

    % eta bounds
    eta_min = zeros(nEta,1);
    eta_max = 50*ones(nEta,1);

    % plant + estimator states
    x_true = zeros(nx,1);
    x_hat_e = zeros(nx_e,1);
    u_prev = zeros(nu,1);

    sigma_y = sqrt(diag(R_meas));

    y_abs = zeros(ny,Tsim);
    u_abs = zeros(nu,Tsim);
    d_hat = zeros(nd,Tsim);
    eps0  = zeros(ny,Tsim);

    % NEW:
    y_true_abs = zeros(ny,Tsim);
    y_meas_abs = zeros(ny,Tsim);

    for k=1:Tsim
        if k>1
            x_true = A*x_true + B*u_prev + E*d_dev(:,k-1);
        end

        y_dev = C*x_true;
        y_abs(:,k) = zs + y_dev;
        y_meas = y_dev + sigma_y.*randn(ny,1);

        % NEW store true/meas ABS
        y_true_abs(:,k) = y_abs(:,k);
        y_meas_abs(:,k) = zs + y_meas;

        innov = y_meas - C_e*x_hat_e;
        x_hat_e = x_hat_e + K_e*innov;
        d_hat(:,k) = x_hat_e(end-nd+1:end);

        Zbar = stackRef(z_ref_dev, k, N);       % ny*N stacked
        y0 = Phi_z_x*x_hat_e;

        f_z  = Phi_z_u' * Wz_big * (y0 - Zbar);

        Ed_prev = [u_prev; zeros((N-1)*nu,1)];
        f_du = -E_du' * Wdu_big * Ed_prev;

        gU = f_z + f_du;
        f = [gU; zeros(nEta,1)];

        % Bounds
        l = [lb_u_stack; eta_min];
        u = [ub_u_stack; eta_max];

        % Constraints:
        % (1) du bounds: du_min <= E_du*U - Ed_prev <= du_max
        A1  = [E_du, zeros(N*nu,nEta)];
        bl1 = du_min_stack + Ed_prev;
        bu1 = du_max_stack + Ed_prev;

        % (2) soft output upper bound: y0 + Phi_z_u*U <= Zmax + eta
        bigM = 1e6;
        A2 = [Phi_z_u, -eye(nZ)];
        bl2 = -bigM*ones(nZ,1);
        bu2 = Z_max_stack - y0;

        Aqp = [A1; A2];
        bl  = [bl1; bl2];
        bu  = [bu1; bu2];

        xinit = zeros(nU+nEta,1);
        [wopt,~] = problem_7(H,f,l,u,Aqp,bl,bu,xinit);

        Uopt = wopt(1:nU);
        eta  = wopt(nU+1:end);

        u_dev_k = Uopt(1:nu);
        u_k_abs = u_s + u_dev_k;

        u_abs(:,k) = u_k_abs;
        eps0(:,k)  = eta(1:ny);           % first-step slack (ny)

        % estimator time update
        x_hat_e = A_e*x_hat_e + B_e*u_dev_k;

        % store prev input in deviation space
        u_prev = u_dev_k;
    end

    sim = struct();
    sim.y_abs      = y_abs;
    sim.y_true_abs = y_true_abs;     % NEW
    sim.y_meas_abs = y_meas_abs;     % NEW
    sim.u_abs      = u_abs;
    sim.d_hat      = d_hat;
    sim.eps        = eps0;
end

function sim = simP10_SDE(t, x0_abs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
                          A,B,E,C, A_e,B_e,C_e, K_e, Phi_x_e,Phi_u_e, Wz,Wu,Wdu, ...
                          u_s,d_s,zs, z_ref_dev, d_dev, ...
                          lb_u_stack,ub_u_stack, du_min_stack,du_max_stack, Z_max_stack, Tsim)
    %#ok<INUSD>
    Ts = t(2)-t(1);

    ny = numel(output_index);
    nu = 2; nd = 2;
    nx_e = size(A_e,1);
    N = numel(lb_u_stack)/nu;

    Phi_z_x = kron(eye(N),C_e) * Phi_x_e;
    Phi_z_u = kron(eye(N),C_e) * Phi_u_e;

    Wz_big  = kron(eye(N),Wz);
    Wu_big  = kron(eye(N),Wu);
    Wdu_big = kron(eye(N),Wdu);

    E_du = makeEdU(N,nu);

    nU = N*nu;
    nZ = N*ny;
    nEta = nZ;

    w_eta = 1e3;
    Weta = w_eta*eye(nEta);

    H_z  = Phi_z_u' * Wz_big  * Phi_z_u;
    H_u  = Wu_big;
    H_du = E_du'   * Wdu_big * E_du;
    H_UU = H_z + H_u + H_du;
    H = blkdiag(H_UU, Weta);

    eta_min = zeros(nEta,1);
    eta_max = 50*ones(nEta,1);

    tf = t(end);
    [~,~,dW3] = ScalarStdWienerProcess(tf, numel(t), 1);
    [~,~,dW4] = ScalarStdWienerProcess(tf, numel(t), 1);

    x_true = zeros(4,Tsim+1);
    x_true(:,1) = x0_abs;

    x_hat_e = zeros(nx_e,1);
    u_prev = zeros(nu,1);

    y_abs = zeros(ny,Tsim);
    u_abs = zeros(nu,Tsim);
    d_hat = zeros(nd,Tsim);
    eps0  = zeros(ny,Tsim);

    % NEW:
    y_true_abs = zeros(ny,Tsim);
    y_meas_abs_store = zeros(ny,Tsim);

    for k=1:Tsim
        y_true_abs4 = FourTankSystemSensor(x_true(:,k), p);
        y_meas_abs  = y_true_abs4(output_index) + sigma_meas_4(output_index).*randn(ny,1);
        y_abs(:,k) = y_meas_abs;

        % NEW store:
        y_true_abs(:,k) = y_true_abs4(output_index);
        y_meas_abs_store(:,k) = y_meas_abs;

        y_meas_dev = y_meas_abs - zs;

        innov = y_meas_dev - C_e*x_hat_e;
        x_hat_e = x_hat_e + K_e*innov;
        d_hat(:,k) = x_hat_e(end-nd+1:end);

        Zbar = stackRef(z_ref_dev, k, N);

        y0  = Phi_z_x*x_hat_e;
        f_z = Phi_z_u' * Wz_big * (y0 - Zbar);

        Ed_prev = [u_prev; zeros((N-1)*nu,1)];
        f_du = -E_du' * Wdu_big * Ed_prev;

        gU = f_z + f_du;
        f = [gU; zeros(nEta,1)];

        l = [lb_u_stack; eta_min];
        u = [ub_u_stack; eta_max];

        A1  = [E_du, zeros(N*nu,nEta)];
        bl1 = du_min_stack + Ed_prev;
        bu1 = du_max_stack + Ed_prev;

        bigM = 1e6;
        A2 = [Phi_z_u, -eye(nZ)];
        bl2 = -bigM*ones(nZ,1);
        bu2 = Z_max_stack - y0;

        Aqp = [A1; A2];
        bl  = [bl1; bl2];
        bu  = [bu1; bu2];

        xinit = zeros(nU+nEta,1);
        [wopt,~] = problem_7(H,f,l,u,Aqp,bl,bu,xinit);

        Uopt = wopt(1:nU);
        eta  = wopt(nU+1:end);

        u_dev_k = Uopt(1:nu);
        u_k_abs = u_s + u_dev_k;

        u_abs(:,k) = u_k_abs;
        eps0(:,k)  = eta(1:ny);

        % estimator time update
        x_hat_e = A_e*x_hat_e + B_e*u_dev_k;
        u_prev  = u_dev_k;

        % plant propagation
        if k<Tsim
            d_abs_k = d_s + d_dev(:,k);
            fplant = FourTankSystemModified(0, x_true(:,k), u_k_abs, d_abs_k, p);
            noise  = [0;0; sigma3*dW3(k); sigma4*dW4(k)];
            x_true(:,k+1) = x_true(:,k) + Ts*fplant + noise;
        end
    end

    sim = struct();
    sim.y_abs      = y_abs;
    sim.y_true_abs = y_true_abs;          % NEW
    sim.y_meas_abs = y_meas_abs_store;    % NEW
    sim.u_abs      = u_abs;
    sim.d_hat      = d_hat;
    sim.eps        = eps0;
end

%% ------------------------------ utilities --------------------------------

function Z = stackRef(z_ref_dev, k, N)
    % z_ref_dev is ny x Tsim, build stacked vector (ny*N x 1)
    [ny,Tsim] = size(z_ref_dev);
    Z = zeros(ny*N,1);
    for j=0:N-1
        idx = min(k+j, Tsim);
        Z(j*ny+(1:ny)) = z_ref_dev(:,idx);
    end
end

function E_du = makeEdU(N,nu)
    E_du = zeros(N*nu, N*nu);
    E_du(1:nu,1:nu) = eye(nu);
    for j=2:N
        rr = (j-1)*nu + (1:nu);
        cc = (j-1)*nu + (1:nu);
        cc0 = (j-2)*nu + (1:nu);
        E_du(rr,cc)  = eye(nu);
        E_du(rr,cc0) = -eye(nu);
    end
end

function addStepMarker(ax, t_step, labelStr)
    if isempty(t_step), return; end
    yl = ylim(ax);
    line(ax, [t_step t_step], yl, 'LineStyle','--','LineWidth',1.2, ...
        'Color',[1 0 0], 'HandleVisibility','off');
    text(ax, t_step+0.01*(ax.XLim(2)-ax.XLim(1)), yl(2)-0.55*(yl(2)-yl(1)), ...
        labelStr, 'Color',[1 0 0], 'FontSize',12, ...
        'Interpreter','none','Clipping','on');
end

function plotCompareOutputsInputs(t, z_ref_abs, s8, s9, s10, u_s, umin, umax, t_step, figName, bigTitle, mkStep)
    figure('Name',figName);
    sgtitle(bigTitle,'Interpreter','none');

    mkIdx = 1:mkStep:numel(t);

    % --- h1
    subplot(3,1,1); hold on; grid on;
    plot(t, z_ref_abs(1,:), 'k--','LineWidth',1.2,'DisplayName','ref h1');

    % true lines
    plot(t, s8.y_true_abs(1,:),  'b','LineWidth',1.4,'DisplayName','P8 true');
    plot(t, s9.y_true_abs(1,:),  'r','LineWidth',1.4,'DisplayName','P9 true');
    plot(t, s10.y_true_abs(1,:), 'g','LineWidth',1.4,'DisplayName','P10 true');

    % measured crosses (sparse)
    plot(t(mkIdx), s8.y_meas_abs(1,mkIdx),  'bx','LineWidth',1.0,'DisplayName','P8 meas');
    plot(t(mkIdx), s9.y_meas_abs(1,mkIdx),  'rx','LineWidth',1.0,'DisplayName','P9 meas');
    plot(t(mkIdx), s10.y_meas_abs(1,mkIdx), 'gx','LineWidth',1.0,'DisplayName','P10 meas');

    addStepMarker(gca, t_step, 'step in d');
    ylabel('$h_1$ [cm]','Interpreter','latex');
    xlim([0 t(end)]);
    legend('Location','best');

    % --- h2
    subplot(3,1,2); hold on; grid on;
    plot(t, z_ref_abs(2,:), 'k--','LineWidth',1.2,'DisplayName','ref h2');

    plot(t, s8.y_true_abs(2,:),  'b','LineWidth',1.4,'DisplayName','P8 true');
    plot(t, s9.y_true_abs(2,:),  'r','LineWidth',1.4,'DisplayName','P9 true');
    plot(t, s10.y_true_abs(2,:), 'g','LineWidth',1.4,'DisplayName','P10 true');

    plot(t(mkIdx), s8.y_meas_abs(2,mkIdx),  'bx','LineWidth',1.0,'DisplayName','P8 meas');
    plot(t(mkIdx), s9.y_meas_abs(2,mkIdx),  'rx','LineWidth',1.0,'DisplayName','P9 meas');
    plot(t(mkIdx), s10.y_meas_abs(2,mkIdx), 'gx','LineWidth',1.0,'DisplayName','P10 meas');

    addStepMarker(gca, t_step, 'step in d');
    ylabel('$h_2$ [cm]','Interpreter','latex');
    xlim([0 t(end)]);
    legend('Location','best');

    % --- inputs
    subplot(3,1,3); hold on; grid on;
    plot(t, s8.u_abs(1,:),  'b','LineWidth',1.2,'DisplayName','P8 u1');
    plot(t, s8.u_abs(2,:),  'b--','LineWidth',1.2,'DisplayName','P8 u2');
    plot(t, s9.u_abs(1,:),  'r','LineWidth',1.2,'DisplayName','P9 u1');
    plot(t, s9.u_abs(2,:),  'r--','LineWidth',1.2,'DisplayName','P9 u2');
    plot(t, s10.u_abs(1,:), 'g','LineWidth',1.2,'DisplayName','P10 u1');
    plot(t, s10.u_abs(2,:), 'g--','LineWidth',1.2,'DisplayName','P10 u2');

    yline(u_s(1),'k--','LineWidth',1.0,'HandleVisibility','off');
    yline(u_s(2),'k--','LineWidth',1.0,'HandleVisibility','off');
    yline(umin(1),'k:','LineWidth',1.0,'DisplayName','u min');
    yline(umax(1),'k:','LineWidth',1.0,'DisplayName','u max');

    addStepMarker(gca, t_step, 'step in d');
    ylabel('$u$ [cm$^3$/s]','Interpreter','latex');
    xlabel('Time [s]');
    xlim([0 t(end)]);
    legend('Location','best');
end

function plotCompareDisturbances(t, d_dev_true, s8, s9, s10, t_step, figName, bigTitle)
    figure('Name',figName);
    sgtitle(bigTitle,'Interpreter','none');

    subplot(2,1,1); hold on; grid on;
    plot(t, d_dev_true(1,:), 'k--','LineWidth',1.2,'DisplayName','true d1');
    plot(t, s8.d_hat(1,:),   'b','LineWidth',1.2,'DisplayName','P8 hat');
    plot(t, s9.d_hat(1,:),   'r','LineWidth',1.2,'DisplayName','P9 hat');
    plot(t, s10.d_hat(1,:),  'g','LineWidth',1.2,'DisplayName','P10 hat');
    addStepMarker(gca, t_step, 'step in d');
    ylabel('$\Delta F_3$ [cm$^3$/s]','Interpreter','latex');
    xlim([0 t(end)]);
    legend('Location','best');

    subplot(2,1,2); hold on; grid on;
    plot(t, d_dev_true(2,:), 'k--','LineWidth',1.2,'DisplayName','true d2');
    plot(t, s8.d_hat(2,:),   'b','LineWidth',1.2,'DisplayName','P8 hat');
    plot(t, s9.d_hat(2,:),   'r','LineWidth',1.2,'DisplayName','P9 hat');
    plot(t, s10.d_hat(2,:),  'g','LineWidth',1.2,'DisplayName','P10 hat');
    addStepMarker(gca, t_step, 'step in d');
    ylabel('$\Delta F_4$ [cm$^3$/s]','Interpreter','latex');
    xlabel('Time [s]');
    xlim([0 t(end)]);
    legend('Location','best');
end

function plotP10Slack(t, s10_lin, s10_sde, t_step, figName, bigTitle)
    figure('Name',figName);
    sgtitle(bigTitle,'Interpreter','none');

    subplot(2,1,1); hold on; grid on;
    plot(t, s10_lin.eps(1,:), 'g','LineWidth',1.4,'DisplayName','P10 eps1 (lin)');
    plot(t, s10_sde.eps(1,:), 'g--','LineWidth',1.4,'DisplayName','P10 eps1 (SDE)');
    addStepMarker(gca, t_step, 'step in d');
    ylabel('$\epsilon_1$','Interpreter','latex');
    xlim([0 t(end)]);
    legend('Location','best');

    subplot(2,1,2); hold on; grid on;
    plot(t, s10_lin.eps(2,:), 'g','LineWidth',1.4,'DisplayName','P10 eps2 (lin)');
    plot(t, s10_sde.eps(2,:), 'g--','LineWidth',1.4,'DisplayName','P10 eps2 (SDE)');
    addStepMarker(gca, t_step, 'step in d');
    ylabel('$\epsilon_2$','Interpreter','latex');
    xlabel('Time [s]');
    xlim([0 t(end)]);
    legend('Location','best');
end

function saveAllFiguresAsPDF(figDir)
    figs = findall(0,'Type','figure');
    for i = 1:numel(figs)
        fig = figs(i);
        name = get(fig,'Name');
        if isempty(name), name = sprintf('Figure_%d', fig.Number); end
        name = regexprep(name,'[^a-zA-Z0-9_\-]','_');
        outFile = fullfile(figDir, [name '.pdf']);

        try
            exportgraphics(fig, outFile, 'ContentType','vector');
        catch
            print(fig, outFile, '-dpdf', '-painters');
        end
    end
end
