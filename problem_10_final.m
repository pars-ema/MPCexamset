%% ========================================================================
% PROBLEM 10 – Input + Soft Output Constrained OFFSET-FREE MPC
% P8-style plots + overlay linear vs SDE plant
%
% Layout style:
%   - Exp-like: outputs figure (h1,h2,u), disturbances figure (d-hat),
%     plus (optional) slack figure.
%   - Legends: Location='best' (like your P8)
%   - Step annotation: red dashed vertical line + text "step in d"
%
% Requires on path:
%   - FourTankSystemModified.m
%   - FourTankSystemSensor.m
%   - ScalarStdWienerProcess.m
%   - problem_7.m
% ========================================================================

clear; clc; close all;

%% ========================================================================
% Global figure style (match Problem 8 robustness)
% ========================================================================
figDir = fullfile('figures','problem_10');
if ~exist(figDir,'dir'), mkdir(figDir); end

set(groot,'defaultFigureColor','w');
set(groot,'defaultFigurePaperPositionMode','auto');
set(groot,'defaultFigurePosition',[100 100 1200 750]);

set(groot,'defaultTextInterpreter','none');
set(groot,'defaultLegendInterpreter','none');
set(groot,'defaultAxesTickLabelInterpreter','none');

set(groot,'defaultAxesFontSize',13);
set(groot,'defaultLegendFontSize',12);

%% ========================================================================
% 0) Physical parameters + steady-state
% ========================================================================
Ts = 10;

a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327;
g  = 981;
rho = 1;
gamma1 = 0.58;
gamma2 = 0.72;

p = [a1;a2;a3;a4;A1;A2;A3;A4;g;gamma1;gamma2;rho];

u_s = [300;300];
d_s = [250;250];

xs_guess = 500*ones(4,1);
opts = optimoptions('fsolve','Display','none');
xs = fsolve(@(x)FourTankSystemModified(0,x,u_s,d_s,p), xs_guess, opts);

h_s_abs_all = FourTankSystemSensor(xs,p);
output_index = [1 2];
zs = h_s_abs_all(output_index);

fprintf("xs:\n"); disp(xs);
fprintf("zs:\n"); disp(zs);

%% ========================================================================
% 1) Linearization + discretization (same as P8/P9)
% ========================================================================
nx = 4; nu = 2; nd = 2; ny = 2;

a1=p(1); a2=p(2); a3=p(3); a4=p(4);
A1=p(5); A2=p(6); A3=p(7); A4=p(8);
g =p(9); gamma1=p(10); gamma2=p(11); rho=p(12);

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

M = [A_c, [B_c E_c];
     zeros(nu+nd, nx+nu+nd)];
Md = expm(M*Ts);

A = Md(1:nx,1:nx);
BEd = Md(1:nx,nx+1:end);
B = BEd(:,1:nu);
E = BEd(:,nu+1:end);
C = C_c;

G = E;

%% ========================================================================
% 2) Noise covariances + augmented model + KF
% ========================================================================
Qw    = diag([25 25]);
Qx    = G * Qw * G';
Qd_RW = 1 * eye(nd);
R_meas = diag([4 4]);

A_e = [A, E;
       zeros(nd,nx), eye(nd)];
B_e = [B;
       zeros(nd,nu)];
C_e = [C, zeros(ny,nd)];
nx_e = nx + nd;

Q_e = blkdiag(Qx, Qd_RW);
K_e = computeKF(A_e,C_e,Q_e,R_meas);

fprintf("KF gain norm: %.4f\n", norm(K_e));

%% ========================================================================
% 3) MPC prediction matrices
% ========================================================================
N = 15;

Wz  = diag([1 1]);
Wu  = 1e-4 * eye(nu);
Wdu = 1e-7 * eye(nu);

[Phi_x_e, Phi_u_e] = buildPredictionMatrices(A_e,B_e,N);

Phi_z_x = kron(eye(N),C_e) * Phi_x_e;
Phi_z_u = kron(eye(N),C_e) * Phi_u_e;

Wz_big  = kron(eye(N),Wz);
Wu_big  = kron(eye(N),Wu);
Wdu_big = kron(eye(N),Wdu);

E_du = zeros(N*nu, N*nu);
E_du(1:nu,1:nu) = eye(nu);
for j = 2:N
    row = (j-1)*nu + (1:nu);
    col = (j-1)*nu + (1:nu);
    col_prev = (j-2)*nu + (1:nu);
    E_du(row,col)      = eye(nu);
    E_du(row,col_prev) = -eye(nu);
end

%% ========================================================================
% 4) Hard input constraints + move constraints
% ========================================================================
u_min_abs = [0;0];
u_max_abs = [600;600];

u_min_dev = u_min_abs - u_s;
u_max_dev = u_max_abs - u_s;

du_max_abs = [300;300];
du_min_abs = -du_max_abs;

u_min_stack  = kron(ones(N,1),u_min_dev);
u_max_stack  = kron(ones(N,1),u_max_dev);

du_min_stack = kron(ones(N,1),du_min_abs);
du_max_stack = kron(ones(N,1),du_max_abs);

%% ========================================================================
% 5) Soft output constraints (absolute bounds)
% ========================================================================
z_min_abs = [0;   0  ];
z_max_abs = [150; 150];

z_min_dev = z_min_abs - zs;
z_max_dev = z_max_abs - zs;

Z_max_stack = kron(ones(N,1), z_max_dev);

%% ========================================================================
% 6) Scenario (like P8 exp B: ref steps + disturbance step)
% ========================================================================
Tsim = 300;
t = (0:Tsim-1) * Ts;

t_change_h1 = [0 600 1600];
values_h1   = [zs(1) zs(1)*1.3 zs(1)*0.8];

t_change_h2 = [0 600 1300];
values_h2   = [zs(2) zs(2)*1.3 zs(2)*1.0];

[h1_ref, h2_ref] = generate_stair_references(t, t_change_h1, values_h1, ...
                                                t_change_h2, values_h2);

z_ref_abs = [h1_ref; h2_ref];
z_ref_dev = z_ref_abs - zs;

d_dev = zeros(nd,Tsim);
k_step = round(2000/Ts);
k_step = max(1, min(Tsim, k_step));
d_step_val = [50; -20];
d_dev(:,k_step:end) = repmat(d_step_val,1,Tsim-k_step+1);

%% ========================================================================
% 7) Linear MPC sim + SDE MPC sim
% ========================================================================
sim_lin = simulateMPC_P10_soft_linear( ...
    A,B,E,C, A_e,B_e,C_e, K_e, ...
    Phi_z_x,Phi_z_u, Wz_big,Wu_big,Wdu_big, E_du, ...
    u_s, u_min_stack,u_max_stack, du_min_stack,du_max_stack, ...
    Z_max_stack, z_ref_dev, zs, d_dev, Tsim, Ts, R_meas);

sigma3 = 2; sigma4 = 2;
sigma_meas_4 = [sqrt(R_meas(1,1)); sqrt(R_meas(2,2)); 0; 0];
x0_abs = xs;

sim_sde = simulateMPC_P10_soft_SDEplant( ...
    t, x0_abs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
    A,B,E,C, A_e,B_e,C_e, K_e, ...
    Phi_z_x,Phi_z_u, Wz_big,Wu_big,Wdu_big, E_du, ...
    u_s, d_s, u_min_stack,u_max_stack, du_min_stack,du_max_stack, ...
    Z_max_stack, z_ref_dev, zs, d_dev, Tsim);

%% ========================================================================
% 8) P8-style plots (same layout) + show soft constraints
% ========================================================================

% ---------- Outputs: h1, h2, u ----------
figure('Name','P10_outputs');
sgtitle('Problem 10: ref steps + disturbance step (linear vs SDE, soft bounds)','Interpreter','none');

subplot(3,1,1); hold on; grid on;
plot(t, z_ref_abs(1,:),          'k--','LineWidth',1.2,'DisplayName','h1 ref');
plot(t, sim_lin.y_abs(1,:),      'b',  'LineWidth',1.4,'DisplayName','h1 linear');
plot(t, sim_sde.y_abs(1,:),      'm:', 'LineWidth',1.4,'DisplayName','h1 SDE meas');
yline(z_min_abs(1),'g:','LineWidth',1.2,'DisplayName','h1 min soft');
yline(z_max_abs(1),'g:','LineWidth',1.2,'DisplayName','h1 max soft');
ylim([60 z_max_abs(1)+5])
addStepMarker(gca, t(k_step), 'step in d');
ylabel('$h_1$ [cm]','Interpreter','latex');
xlim([0 t(end)]);
legend("Location","best");

subplot(3,1,2); hold on; grid on;
plot(t, z_ref_abs(2,:),          'k--','LineWidth',1.2,'DisplayName','h2 ref');
plot(t, sim_lin.y_abs(2,:),      'b',  'LineWidth',1.4,'DisplayName','h2 linear');
plot(t, sim_sde.y_abs(2,:),      'm:', 'LineWidth',1.4,'DisplayName','h2 SDE meas');
yline(z_min_abs(2),'g:','LineWidth',1.2,'DisplayName','h2 min soft');
yline(z_max_abs(2),'g:','LineWidth',1.2,'DisplayName','h2 max soft');
ylim([100 160])
addStepMarker(gca, t(k_step), 'step in d');
ylabel('$h_2$ [cm]','Interpreter','latex');
xlim([0 t(end)]);
legend("Location","best");

subplot(3,1,3); hold on; grid on;
plot(t, sim_lin.u_abs(1,:), 'b','LineWidth',1.4,'DisplayName','u1 linear');
plot(t, sim_lin.u_abs(2,:), 'r','LineWidth',1.4,'DisplayName','u2 linear');
plot(t, sim_sde.u_abs(1,:), 'm:','LineWidth',1.2,'DisplayName','u1 SDE');
plot(t, sim_sde.u_abs(2,:), 'm-.','LineWidth',1.2,'DisplayName','u2 SDE');
yline(u_s(1),'k--','LineWidth',1.0,'HandleVisibility','off');
yline(u_s(2),'k--','LineWidth',1.0,'HandleVisibility','off');
yline(u_min_abs(1),'k:','LineWidth',1.0,'DisplayName','u min');
yline(u_max_abs(1),'k:','LineWidth',1.0,'DisplayName','u max');
addStepMarker(gca, t(k_step), 'step in d');
ylabel('$u$ [cm$^3$/s]','Interpreter','latex');
xlabel('Time [s]');
xlim([0 t(end)]);
legend("Location","best");

% ---------- Disturbance estimates ----------
figure('Name','P10_disturbances');
sgtitle('Problem 10: disturbance estimates (linear vs SDE)','Interpreter','none');

subplot(2,1,1); hold on; grid on;
plot(t, d_dev(1,:),           'k--','LineWidth',1.2,'DisplayName','d1 true');
plot(t, sim_lin.d_hat(1,:),   'b',  'LineWidth',1.4,'DisplayName','d1 hat linear');
plot(t, sim_sde.d_hat(1,:),   'm:', 'LineWidth',1.2,'DisplayName','d1 hat SDE');
addStepMarker(gca, t(k_step), 'step in d');
ylabel('$\Delta F_3$ [cm$^3$/s]','Interpreter','latex');
xlim([0 t(end)]);
legend("Location","best");

subplot(2,1,2); hold on; grid on;
plot(t, d_dev(2,:),           'k--','LineWidth',1.2,'DisplayName','d2 true');
plot(t, sim_lin.d_hat(2,:),   'b',  'LineWidth',1.4,'DisplayName','d2 hat linear');
plot(t, sim_sde.d_hat(2,:),   'm:', 'LineWidth',1.2,'DisplayName','d2 hat SDE');
addStepMarker(gca, t(k_step), 'step in d');
ylabel('$\Delta F_4$ [cm$^3$/s]','Interpreter','latex');
xlabel('Time [s]');
xlim([0 t(end)]);
legend("Location","best");

% ---------- Slack (optional but very useful) ----------
figure('Name','P10_slacks');
sgtitle('Problem 10: first-step slack (soft constraint activity)','Interpreter','none');

subplot(2,1,1); hold on; grid on;
plot(t, sim_lin.eps(1,:), 'b','LineWidth',1.4,'DisplayName','eps1 linear');
plot(t, sim_sde.eps(1,:), 'm:','LineWidth',1.2,'DisplayName','eps1 SDE');
addStepMarker(gca, t(k_step), 'step in d');
ylabel('$\epsilon_1$ [cm]','Interpreter','latex');
xlim([0 t(end)]);
legend("Location","best");

subplot(2,1,2); hold on; grid on;
plot(t, sim_lin.eps(2,:), 'b','LineWidth',1.4,'DisplayName','eps2 linear');
plot(t, sim_sde.eps(2,:), 'm:','LineWidth',1.2,'DisplayName','eps2 SDE');
addStepMarker(gca, t(k_step), 'step in d');
ylabel('$\epsilon_2$ [cm]','Interpreter','latex');
xlabel('Time [s]');
xlim([0 t(end)]);
legend("Location","best");

%% ========================================================================
% 9) Save all figures as vector PDFs
% ========================================================================
saveAllFiguresAsPDF(figDir);
fprintf('Saved figures to: %s\n', figDir);

%% ========================================================================
% ============================ FUNCTIONS ==================================
% ========================================================================

function K = computeKF(A,C,Q,R)
    [P,~,~] = dare(A',C',Q,R);
    P = 0.5*(P+P');
    Re = C*P*C' + R;
    K  = P*C'/Re;
end

function [Phi_x, Gamma_x] = buildPredictionMatrices(A,B,N)
    [nx,nu] = size(B);
    Phi_x   = zeros(N*nx, nx);
    Gamma_x = zeros(N*nx, N*nu);

    A_power = eye(nx);
    for i = 1:N
        A_power = A * A_power;
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

function sim = simulateMPC_P10_soft_linear( ...
    A,B,E,C, A_e,B_e,C_e, K_e, ...
    Phi_z_x,Phi_z_u, Wz_big,Wu_big,Wdu_big, E_du, ...
    u_s, u_min_stack,u_max_stack, du_min_stack,du_max_stack, ...
    Z_max_stack, z_ref_dev, zs, d_dev, Tsim, Ts, R_meas)

    [nx,nu] = size(B);
    nd = size(E,2);
    ny = size(C,1);
    nx_e = nx + nd;

    N  = size(Phi_z_u,1) / ny;
    nU = N*nu;
    nZ = N*ny;
    nEta = nZ;

    % Slack penalty
    w_eta = 1e3;
    Weta_big = w_eta * eye(nEta);

    % Hessian blocks
    H_z  = Phi_z_u' * Wz_big  * Phi_z_u;
    H_u  = Wu_big;
    H_du = E_du'   * Wdu_big * E_du;
    H_UU = H_z + H_u + H_du;

    % eta bounds
    eta_max = 50;
    eta_min_vec = zeros(nEta,1);
    eta_max_vec = eta_max * ones(nEta,1);

    % plant noise (linear)
    Qplant = 1e-4 * eye(nx);
    sigma_meas = sqrt(diag(R_meas));

    % init
    x_true = zeros(nx,1);
    x_hat_e = zeros(nx_e,1);
    u_prev = zeros(nu,1);

    y_abs = zeros(ny,Tsim);
    u_abs = zeros(nu,Tsim);
    d_hat = zeros(nd,Tsim);
    eps0  = zeros(ny,Tsim);

    nu_total = nU + nEta;

    for k = 1:Tsim

        if k > 1
            w = mvnrnd(zeros(nx,1),Qplant)'; %#ok<MVNRND>
            x_true = A*x_true + B*u_prev + E*d_dev(:,k-1) + w;
        end

        y_dev = C*x_true;
        y_abs(:,k) = zs + y_dev;

        y_meas = y_dev + sigma_meas .* randn(ny,1);

        % KF update
        innov  = y_meas - C_e*x_hat_e;
        x_hat_e = x_hat_e + K_e*innov;
        d_hat(:,k) = x_hat_e(nx+1:end);

        % horizon ref
        if k+N-1 <= Tsim
            Zbar = z_ref_dev(:,k:k+N-1);
        else
            last = Tsim-k+1;
            Zbar = [z_ref_dev(:,k:end), repmat(z_ref_dev(:,end),1,N-last)];
        end
        Zbar = Zbar(:);

        y0 = Phi_z_x*x_hat_e;

        f_z = Phi_z_u' * Wz_big * (y0 - Zbar);

        Ed_prev = [u_prev; zeros((N-1)*nu,1)];
        f_du = -E_du' * Wdu_big * Ed_prev;

        g_U = f_z + f_du;

        H = blkdiag(H_UU, Weta_big);
        f = [g_U; zeros(nEta,1)];

        l = [u_min_stack; eta_min_vec];
        u = [u_max_stack; eta_max_vec];

        % move constraints
        A1  = [E_du, zeros(N*nu,nEta)];
        bl1 = du_min_stack + Ed_prev;
        bu1 = du_max_stack + Ed_prev;

        % soft upper bound: y0 + Phi_u U <= Zmax + eta
        bigM = 1e6;
        b_soft = Z_max_stack - y0;
        A2  = [Phi_z_u, -eye(nZ)];
        bl2 = -bigM*ones(nZ,1);
        bu2 = b_soft;

        Aqp = [A1; A2];
        bl  = [bl1; bl2];
        bu  = [bu1; bu2];

        xinit = zeros(nu_total,1);
        [w_opt,~] = problem_7(H,f,l,u,Aqp,bl,bu,xinit);

        Uopt = w_opt(1:nU);
        etaopt = w_opt(nU+1:end);

        u_k_dev = Uopt(1:nu);
        u_k_abs = u_s + u_k_dev;

        u_abs(:,k) = u_k_abs;
        eps0(:,k)  = etaopt(1:ny);

        % KF time update
        x_hat_e = A_e*x_hat_e + B_e*u_k_dev;
        u_prev  = u_k_dev;
    end

    sim = struct();
    sim.y_abs = y_abs;
    sim.u_abs = u_abs;
    sim.d_hat = d_hat;
    sim.eps   = eps0;
end

function sim = simulateMPC_P10_soft_SDEplant( ...
    t, x0_abs, xs, p, output_index, sigma_meas_4, sigma3, sigma4, ...
    A,B,E,C, A_e,B_e,C_e, K_e, ...
    Phi_z_x,Phi_z_u, Wz_big,Wu_big,Wdu_big, E_du, ...
    u_s, d_s, u_min_stack,u_max_stack, du_min_stack,du_max_stack, ...
    Z_max_stack, z_ref_dev, zs, d_dev, Tsim)

    %#ok<INUSD>
    Ts = t(2)-t(1);

    ny = numel(output_index);
    [~,nu] = size(B);
    nd = size(E,2);
    nx_e = size(A_e,1);

    N  = size(Phi_z_u,1) / ny;
    nU = N*nu;
    nZ = N*ny;
    nEta = nZ;

    w_eta = 1e3;
    Weta_big = w_eta * eye(nEta);

    H_z  = Phi_z_u' * Wz_big  * Phi_z_u;
    H_u  = Wu_big;
    H_du = E_du'   * Wdu_big * E_du;
    H_UU = H_z + H_u + H_du;

    eta_max = 50;
    eta_min_vec = zeros(nEta,1);
    eta_max_vec = eta_max * ones(nEta,1);

    tf = t(end);
    [~,~,dW3] = ScalarStdWienerProcess(tf, numel(t), 1);
    [~,~,dW4] = ScalarStdWienerProcess(tf, numel(t), 1);

    x_true_abs = zeros(4,Tsim+1);
    x_true_abs(:,1) = x0_abs;

    x_hat_e = zeros(nx_e,1);
    u_prev = zeros(nu,1);

    y_abs = zeros(ny,Tsim);
    u_abs = zeros(nu,Tsim);
    d_hat = zeros(nd,Tsim);
    eps0  = zeros(ny,Tsim);

    nu_total = nU + nEta;

    for k = 1:Tsim

        y_true_abs4 = FourTankSystemSensor(x_true_abs(:,k), p);
        y_meas_abs  = y_true_abs4(output_index) + sigma_meas_4(output_index).*randn(ny,1);
        y_meas_dev  = y_meas_abs - zs;

        y_abs(:,k) = y_meas_abs;

        innov  = y_meas_dev - C_e*x_hat_e;
        x_hat_e = x_hat_e + K_e*innov;
        d_hat(:,k) = x_hat_e(end-nd+1:end);

        if k+N-1 <= Tsim
            Zbar = z_ref_dev(:,k:k+N-1);
        else
            last = Tsim-k+1;
            Zbar = [z_ref_dev(:,k:end), repmat(z_ref_dev(:,end),1,N-last)];
        end
        Zbar = Zbar(:);

        y0 = Phi_z_x*x_hat_e;
        f_z = Phi_z_u' * Wz_big * (y0 - Zbar);

        Ed_prev = [u_prev; zeros((N-1)*nu,1)];
        f_du = -E_du' * Wdu_big * Ed_prev;

        g_U = f_z + f_du;

        H = blkdiag(H_UU, Weta_big);
        f = [g_U; zeros(nEta,1)];

        l = [u_min_stack; eta_min_vec];
        u = [u_max_stack; eta_max_vec];

        A1  = [E_du, zeros(N*nu,nEta)];
        bl1 = du_min_stack + Ed_prev;
        bu1 = du_max_stack + Ed_prev;

        bigM = 1e6;
        b_soft = Z_max_stack - y0;
        A2  = [Phi_z_u, -eye(nZ)];
        bl2 = -bigM*ones(nZ,1);
        bu2 = b_soft;

        Aqp = [A1; A2];
        bl  = [bl1; bl2];
        bu  = [bu1; bu2];

        xinit = zeros(nu_total,1);
        [w_opt,~] = problem_7(H,f,l,u,Aqp,bl,bu,xinit);

        Uopt = w_opt(1:nU);
        etaopt = w_opt(nU+1:end);

        u_k_dev = Uopt(1:nu);
        u_k_abs = u_s + u_k_dev;

        u_abs(:,k) = u_k_abs;
        eps0(:,k)  = etaopt(1:ny);

        x_hat_e = A_e*x_hat_e + B_e*u_k_dev;
        u_prev  = u_k_dev;

        if k < Tsim
            d_abs_k = d_s + d_dev(:,k);
            fplant = FourTankSystemModified(0, x_true_abs(:,k), u_k_abs, d_abs_k, p);

            noise = [0; 0; sigma3*dW3(k); sigma4*dW4(k)];
            x_true_abs(:,k+1) = x_true_abs(:,k) + Ts*fplant + noise;
        end
    end

    sim = struct();
    sim.y_abs = y_abs;
    sim.u_abs = u_abs;
    sim.d_hat = d_hat;
    sim.eps   = eps0;
end

function addStepMarker(ax, xStep, labelStr)
% matches your P8 manual line + text method
    yl = ylim(ax);
    line(ax, [xStep xStep], yl, 'LineStyle','--','LineWidth',1.2, ...
        'Color',[1 0 0], 'HandleVisibility','off');
    text(ax, xStep+0.01*(ax.XLim(2)-ax.XLim(1)), yl(2)-0.5*(yl(2)-yl(1)), ...
        labelStr, 'Color',[1 0 0], 'FontSize',12, ...
        'Interpreter','none','Clipping','on');
end

function [h1_ref, h2_ref] = generate_stair_references(t, t_change_h1, values_h1, t_change_h2, values_h2)
    h1_ref = zeros(1,length(t));
    h2_ref = zeros(1,length(t));

    for k = 1:length(t)
        idx = find(t(k) >= t_change_h1, 1, 'last');
        if isempty(idx), idx = 1; end
        idx = min(idx, numel(values_h1));
        h1_ref(k) = values_h1(idx);
    end

    for k = 1:length(t)
        idx = find(t(k) >= t_change_h2, 1, 'last');
        if isempty(idx), idx = 1; end
        idx = min(idx, numel(values_h2));
        h2_ref(k) = values_h2(idx);
    end
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
