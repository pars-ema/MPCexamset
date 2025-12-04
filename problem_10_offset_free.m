%% ========================================================================
% PROBLEM 10 – Input + Soft Output Constrained OFFSET-FREE MPC
% Using SAME disturbance augmentation + KF idea as Problems 8 & 9:
%
%   x_{k+1} = A x_k + B u_dev_k + E d_k
%   d_{k+1} = d_k                  (random walk)
%   y_k     = C x_k
%
% Augmented state for KF + MPC: x_e = [x; d]
%
% New vs Problem 9:
%   - SAME hard input constraints as before
%   - PLUS soft output constraints on h1, h2 via slack variables
%
% Decision variables in QP:
%   w = [U_seq; eps_seq]
%     U_seq : stacked input deviations (N*nu x 1)
%     eps   : stacked slacks for outputs (N*ny x 1)
%
% QP solved with problem_7.m:
%   [w_opt, info] = problem_7(H,f,l,u,Aqp,bl,bu,xinit)
% ========================================================================

clear; clc; close all;

%% ========================================================================
% 0) Physical parameters + steady-state (same as P6/P8/P9)
% ========================================================================

Ts = 10;                % Sampling time [s]

a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;        % [cm^2]
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327;% [cm^2]
g  = 981;                                                   % [cm/s^2]
rho = 1;                                                    % [g/cm^3]
gamma1 = 0.58;
gamma2 = 0.72;

p = [a1;a2;a3;a4;A1;A2;A3;A4;g;gamma1;gamma2;rho];

% Operating point (same as P6/8/9)
u_s = [300;300];          % [cm^3/s] (F1, F2)
d_s = [250;250];          % [cm^3/s] (F3, F4)

% Steady-state masses xs (for nonlinear model / physical linearization)
xs_guess = 500*ones(4,1);
opts = optimoptions('fsolve','Display','none');
xs = fsolve(@(x)FourTankSystemModified(0,x,u_s,d_s,p), xs_guess, opts);

fprintf("Steady-state xs (masses):\n");
disp(xs);

% Steady-state heights (all 4 tanks)
h_s_abs_all = FourTankSystemSensor(xs,p);   % [4x1] heights in cm
output_index = [1 2];                       % we measure tanks 1 and 2
zs = h_s_abs_all(output_index);             % steady-state measured heights [2x1]

fprintf("Steady-state zs = [h1_s; h2_s]:\n");
disp(zs);

%% ========================================================================
% 1) Continuous linearization & discretization (physical model, deviation)
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

% Continuous-time A matrix (mass states, deviation form)
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

% Continuous-time disturbance input matrix E_c (F3,F4 → tanks 3 & 4)
E_c = [
    0, 0;
    0, 0;
    1, 0;
    0, 1
];

% Output matrix: masses → h1,h2
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
% 2) Noise covariances & augmented model (OFFSET-FREE)
% ========================================================================

% Process noise on disturbance channels w_k (2D: affects tanks 3 & 4)
Qw = diag([25 25]);         % same as in KF problems

% Equivalent process noise on states: Qx = G Qw G'
Qx = G * Qw * G';

% Random-walk disturbance model covariance
Qd_RW = 1 * eye(nd);

% Measurement noise covariance R (heights h1,h2)
R_meas = diag([4 4]);       % (std = 2 cm)

% Augmented model x_e = [x; d]:
%   x_{k+1} = A x_k + B u_dev_k + E d_k
%   d_{k+1} = d_k
%   y_k     = C x_k
A_e = [A, E;
       zeros(nd,nx), eye(nd)];
B_e = [B;
       zeros(nd,nu)];
C_e = [C, zeros(ny,nd)];

nx_e = nx + nd;

% Augmented process covariance
Q_e = blkdiag(Qx, Qd_RW);

% Static augmented Kalman filter (stationary K_e)
K_e = computeKF(A_e,C_e,Q_e,R_meas);
fprintf("Static augmented KF gain norm: %.4f\n", norm(K_e));

%% ========================================================================
% 3) MPC prediction matrices (augmented model, deviation form)
% ========================================================================

N = 15;    % prediction horizon

% Weights
Wz  = diag([1 1]);      % output tracking
Wu  = 1e-4 * eye(nu);   % input magnitude
Wdu = 1e-7 * eye(nu);   % input moves

[Phi_x_e, Phi_u_e] = buildPredictionMatrices(A_e,B_e,N);

Phi_z_x = kron(eye(N),C_e) * Phi_x_e;   % (N*ny x nx_e)
Phi_z_u = kron(eye(N),C_e) * Phi_u_e;   % (N*ny x N*nu)

Wz_big  = kron(eye(N),Wz);
Wu_big  = kron(eye(N),Wu);
Wdu_big = kron(eye(N),Wdu);

% Δu operator
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
% 4) Input constraints (absolute → deviation)
% ========================================================================

u_min_abs = [0;0];
u_max_abs = [450;450];

u_min_dev = u_min_abs - u_s;
u_max_dev = u_max_abs - u_s;

% move constraints Δu_min <= Δu <= Δu_max
du_max_abs = [300;300];
du_min_abs = -du_max_abs;

du_max_dev = du_max_abs;
du_min_dev = du_min_abs;

u_min_stack  = kron(ones(N,1),u_min_dev);
u_max_stack  = kron(ones(N,1),u_max_dev);

du_min_stack = kron(ones(N,1),du_min_dev);
du_max_stack = kron(ones(N,1),du_max_dev);

%% ========================================================================
% 5) Soft OUTPUT constraints (absolute → deviation)
% ========================================================================

% Soft bounds on absolute heights (tune as you like)
z_min_abs = [0;   0  ];     % lower soft bounds [cm]
z_max_abs = [150; 150];     % upper soft bounds [cm] ~ tank capacity

% Convert to deviation coordinates: z_dev = z_abs - zs
z_min_dev = z_min_abs - zs;
z_max_dev = z_max_abs - zs;

Z_min_stack = kron(ones(N,1), z_min_dev);   %#ok<NASGU> (not used explicitly)
Z_max_stack = kron(ones(N,1), z_max_dev);

%% ========================================================================
% 6) Simulation setup (reference + disturbance step)
% ========================================================================

Tsim = 300;
t = (0:Tsim-1) * Ts;

% Reference staircases (same style as you had)
t_change_h1 = [0 600 1600];
values_h1   = [zs(1) zs(1)*1.3 zs(1)*0.8];

t_change_h2 = [0 600 1300];
values_h2   = [zs(2) zs(2)*1.3 zs(2)*1.0];

[h1_ref, h2_ref] = generate_stair_references( ...
                        t, t_change_h1, values_h1, ...
                           t_change_h2, values_h2);

z_ref_abs = [h1_ref; h2_ref];    % absolute refs
z_ref     = z_ref_abs - zs;      % deviation refs

% Disturbance scenario (deviation from d_s):
%   0 until t ≈ 2000 s, then step in F3
d_dev = zeros(nd,Tsim);
k_step = round(2000 / Ts);
if k_step < 1, k_step = 1; end
d_step_val = [50; -20];   % [ΔF3; ΔF4] in cm^3/s, quite visible
d_dev(:,k_step:end) = repmat(d_step_val,1,Tsim-k_step+1);

%% ========================================================================
% 7) Closed-loop simulation (OFFSET-FREE MPC + static augmented KF)
% ========================================================================

sim = simulateMPC_Problem10_offsetfree( ...
    A,B,E,C, ...
    A_e,B_e,C_e, ...
    K_e, ...
    Phi_z_x,Phi_z_u, ...
    Wz_big,Wu_big,Wdu_big, ...
    E_du, ...
    u_s, ...
    u_min_stack,u_max_stack, ...
    du_min_stack,du_max_stack, ...
    Z_max_stack, ...
    z_ref_abs,z_ref, ...
    z_min_abs,z_max_abs, ...
    zs, ...
    d_dev, ...
    Tsim,Ts,R_meas);

%% ========================================================================
% 8) Plots – outputs, inputs, constraints, disturbances
% ========================================================================

plotResultsProblem10_offsetfree(t,sim, ...
    z_ref_abs,z_ref, ...
    u_min_abs,u_max_abs, ...
    du_min_abs,du_max_abs, ...
    z_min_abs,z_max_abs, ...
    k_step, d_step_val);

%% ========================================================================
% ======================== LOCAL FUNCTIONS ================================
% ========================================================================

function K = computeKF(A,C,Q,R)
% Stationary (steady-state) Kalman filter gain from DARE
    [P,~,~] = dare(A',C',Q,R);
    P = 0.5*(P+P');               % enforce symmetry
    Re = C*P*C' + R;
    K  = P*C'/Re;
end

% -------------------------------------------------------------------------
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
function sim = simulateMPC_Problem10_offsetfree( ...
    A,B,E,C, ...
    A_e,B_e,C_e, ...
    K_e, ...
    Phi_z_x,Phi_z_u, ...
    Wz_big,Wu_big,Wdu_big, ...
    E_du, ...
    u_s, ...
    u_min_stack,u_max_stack, ...
    du_min_stack,du_max_stack, ...
    Z_max_stack, ...
    z_ref_abs,z_ref, ...
    z_min_abs,z_max_abs, ...
    zs, ...
    d_dev, ...
    Tsim,Ts,R_meas)

    [nx,nu] = size(B);
    nd = size(E,2);
    ny = size(C,1);
    nx_e = nx + nd;

    N  = size(Phi_z_u,1) / ny;   % prediction horizon
    nU = N*nu;
    nZ = N*ny;
    nEta = nZ;                   % one slack per predicted output

    % Slack penalty (softness) – big = stricter
    w_eta = 1e3;
    Weta_big = w_eta * eye(nEta);

    % Precompute constant Hessian part for U
    H_z  = Phi_z_u' * Wz_big  * Phi_z_u;
    H_u  = Wu_big;
    H_du = E_du'   * Wdu_big * E_du;
    H_UU = H_z + H_u + H_du;

    % Box bounds for η (slack)
    eta_max = 50;                    % [cm] extra allowed above z_max_dev
    eta_min_vec = zeros(nEta,1);
    eta_max_vec = eta_max * ones(nEta,1);

    % Plant noise (linear)
    Qplant = 1e-4 * eye(nx);
    sigma_meas = sqrt(diag(R_meas));

    % Deviation bounds for outputs
    z_min_dev = z_min_abs - zs;
    z_max_dev = z_max_abs - zs;
    Z_max_stack_dev = Z_max_stack; %#ok<NASGU> (name clarity)

    % Initialise
    x_true = zeros(nx,1);            % deviation state
    d_true = zeros(nd,1);            % deviation disturbance
    x_hat_e = zeros(nx_e,1);         % augmented estimate [x_hat; d_hat]
    u_prev = zeros(nu,1);

    % Logs
    z_dev = zeros(ny,Tsim);
    z_abs = zeros(ny,Tsim);
    u_dev = zeros(nu,Tsim);
    u_abs = zeros(nu,Tsim);
    d_hat = zeros(nd,Tsim);
    eps0  = zeros(ny,Tsim);          % first-step slack

    nu_total = nU + nEta;

    for k = 1:Tsim

        % Disturbance for this step (deviation)
        d_true = d_dev(:,k);

        % ----- True plant update -----
        if k > 1
            w = mvnrnd(zeros(nx,1),Qplant)';
            x_true = A*x_true + B*u_prev + E*d_dev(:,k-1) + w;
        end

        z_dev(:,k) = C*x_true;
        z_abs(:,k) = zs + z_dev(:,k);

        % Measurement with noise (deviation)
        v = sigma_meas .* randn(ny,1);
        y_meas = z_dev(:,k) + v;

        % ----- Static augmented Kalman filter (offset-free) -----
        % Augmented measurement model: y = C_e x_e, with C_e = [C  0]
        y_pred = C_e*x_hat_e;
        innov  = y_meas - y_pred;
        x_hat_e = x_hat_e + K_e*innov;

        % Store estimated disturbances
        d_hat(:,k) = x_hat_e(nx+1:end);

        % ----- Build horizon references (deviation) -----
        if k+N-1 <= Tsim
            zbar = z_ref(:,k:k+N-1);
            ubar = zeros(nu,N);          % nominal input deviations = 0
        else
            last = Tsim-k+1;
            zbar = [z_ref(:,k:end), repmat(z_ref(:,end),1,N-last)];
            ubar = zeros(nu,N);
        end

        Zbar = zbar(:);
        Ubar = ubar(:);

        % ----- Cost gradient for U -----
        y0 = Phi_z_x*x_hat_e;          % predicted output dev offset term

        f_z  = Phi_z_u' * Wz_big  * (y0 - Zbar);
        f_u  = -Wu_big * Ubar;

        Ed_prev = [u_prev; zeros((N-1)*nu,1)];
        f_du = -E_du' * Wdu_big * Ed_prev;

        g_U = f_z + f_u + f_du;

        % Full decision vector w = [U; eta]
        H = blkdiag(H_UU, Weta_big);
        f = [g_U; zeros(nEta,1)];

        % ----- Box constraints: U and eta -----
        l = [u_min_stack; eta_min_vec];
        u = [u_max_stack; eta_max_vec];

        % ----- Move constraints: ΔUmin <= E_du U - Ed_prev <= ΔUmax -----
        A1  = [E_du, zeros(N*nu,nEta)];
        bl1 = du_min_stack + Ed_prev;
        bu1 = du_max_stack + Ed_prev;

        % ----- Soft output upper bounds: z_pred_dev <= z_max_dev + eta -----
        % z_pred_dev = y0 + Phi_z_u U
        %   y0 + Phi_z_u U <= Z_max_stack + eta
        % => Phi_z_u U - I*eta <= Z_max_stack - y0
        Zmax_stack = Z_max_stack;          % (N*ny x 1) dev-bounds-stack
        b_soft = Zmax_stack - y0;

        bigM = 1e6;
        A2  = [Phi_z_u, -eye(nZ)];
        bl2 = -bigM * ones(nZ,1);         % effectively one-sided
        bu2 = b_soft;

        % Stack inequalities
        Aqp   = [A1; A2];
        bl_ineq = [bl1; bl2];
        bu_ineq = [bu1; bu2];

        % ----- Solve QP with problem_7 -----
        xinit = zeros(nu_total,1);
        [w_opt,~] = problem_7(H,f,l,u,Aqp,bl_ineq,bu_ineq,xinit);

        Uopt   = w_opt(1:nU);
        etaopt = w_opt(nU+1:end);         %#ok<NASGU> (full seq, only first shown)

        u_k = Uopt(1:nu);

        u_dev(:,k) = u_k;
        u_abs(:,k) = u_s + u_k;

        % First-step slack for each output (for plotting)
        eps0(:,k) = etaopt(1:ny);

        % Time update of KF
        x_hat_e = A_e*x_hat_e + B_e*u_k;

        u_prev = u_k;
    end

    % Pack results
    sim.z_dev = z_dev;
    sim.z_abs = z_abs;
    sim.u_dev = u_dev;
    sim.u_abs = u_abs;
    sim.d_hat = d_hat;
    sim.d_true = d_dev;
    sim.eps = eps0;
end

% -------------------------------------------------------------------------
function plotResultsProblem10_offsetfree(t,sim, ...
    z_ref_abs,z_ref, ...
    u_min_abs,u_max_abs, ...
    du_min_abs,du_max_abs, ...
    z_min_abs,z_max_abs, ...
    k_step, d_step_val)

    z_abs = sim.z_abs;
    z_dev = sim.z_dev;
    u_abs = sim.u_abs;
    u_dev = sim.u_dev;
    d_hat = sim.d_hat;
    d_true = sim.d_true;
    eps0  = sim.eps;

    %% Absolute outputs with soft bounds
    figure('Name','P10 – Absolute Outputs (Offset-free, Soft Constraints)');
    subplot(2,1,1); hold on; grid on;
    plot(t, z_abs(1,:), 'b','LineWidth',1.5);
    plot(t, z_ref_abs(1,:), 'k--','LineWidth',1.2);
    yline(z_min_abs(1),'g:','LineWidth',1.2);
    yline(z_max_abs(1),'g:','LineWidth',1.2);
    xline(t(k_step),'r--','d-step','LabelVerticalAlignment','bottom');
    ylabel('h_1 [cm]');
    legend('h_1','h_{1,ref}','h_{1,min}^{soft}','h_{1,max}^{soft}','Location','best');
    title('Tank 1 – Absolute Height with Soft Bounds');

    subplot(2,1,2); hold on; grid on;
    plot(t, z_abs(2,:), 'b','LineWidth',1.5);
    plot(t, z_ref_abs(2,:), 'k--','LineWidth',1.2);
    yline(z_min_abs(2),'g:','LineWidth',1.2);
    yline(z_max_abs(2),'g:','LineWidth',1.2);
    xline(t(k_step),'r--','d-step','LabelVerticalAlignment','bottom');
    ylabel('h_2 [cm]'); xlabel('Time [s]');
    legend('h_2','h_{2,ref}','h_{2,min}^{soft}','h_{2,max}^{soft}','Location','best');
    title('Tank 2 – Absolute Height with Soft Bounds');

    %% Deviation outputs
    figure('Name','P10 – Deviation Outputs (Offset-free)');
    subplot(2,1,1); hold on; grid on;
    plot(t, z_dev(1,:), 'b','LineWidth',1.5);
    plot(t, z_ref(1,:), 'k--','LineWidth',1.2);
    xline(t(k_step),'r--','d-step');
    ylabel('z_1 dev [cm]');
    legend('z_1','z_{1,ref}','Location','best');
    title('Output 1 (deviation)');

    subplot(2,1,2); hold on; grid on;
    plot(t, z_dev(2,:), 'b','LineWidth',1.5);
    plot(t, z_ref(2,:), 'k--','LineWidth',1.2);
    xline(t(k_step),'r--','d-step');
    ylabel('z_2 dev [cm]'); xlabel('Time [s]');
    legend('z_2','z_{2,ref}','Location','best');
    title('Output 2 (deviation)');

    %% Absolute inputs + hard bounds
    figure('Name','P10 – Absolute Inputs (Offset-free)');
    hold on; grid on;
    plot(t, u_abs(1,:), 'r','LineWidth',1.5);
    plot(t, u_abs(2,:), 'm--','LineWidth',1.5);

    yline(u_min_abs(1),'b--','LineWidth',1.2);
    yline(u_max_abs(1),'b--','LineWidth',1.2);
    yline(u_min_abs(2),'g--','LineWidth',1.2);
    yline(u_max_abs(2),'g--','LineWidth',1.2);
    xline(t(k_step),'r--','d-step');

    xlabel('Time [s]');
    ylabel('Flow [cm^3/s]');
    legend('u_1','u_2', ...
           'u_{1,min}','u_{1,max}', ...
           'u_{2,min}','u_{2,max}','Location','best');
    title('Absolute Pump Flows with Input Constraints');

    %% Deviation inputs
    figure('Name','P10 – Deviation Inputs (Offset-free)');
    hold on; grid on;
    stairs(t, u_dev(1,:), 'r','LineWidth',1.5);
    stairs(t, u_dev(2,:), 'm--','LineWidth',1.5);
    xline(t(k_step),'r--','d-step');
    xlabel('Time [s]');
    ylabel('\Delta u [cm^3/s]');
    legend('\Delta u_1','\Delta u_2','Location','best');
    title('Deviation Inputs');

    %% Input moves + move constraints
    du1 = [u_dev(1,1), diff(u_dev(1,:))];
    du2 = [u_dev(2,1), diff(u_dev(2,:))];

    figure('Name','P10 – Input Moves (Offset-free)');
    hold on; grid on;
    stairs(t, du1, 'r','LineWidth',1.5);
    stairs(t, du2, 'm--','LineWidth',1.5);

    yline(du_max_abs(1),'b--');
    yline(du_min_abs(1),'b--');
    yline(du_max_abs(2),'g--');
    yline(du_min_abs(2),'g--');
    xline(t(k_step),'r--','d-step');

    xlabel('Time [s]');
    ylabel('\Delta u [cm^3/s]');
    legend('\Delta u_1','\Delta u_2', ...
           '\Delta u_{1,max}','\Delta u_{1,min}', ...
           '\Delta u_{2,max}','\Delta u_{2,min}','Location','best');
    title('Input Moves and Move Constraints');

    %% Disturbances true vs estimated (offset-free behaviour)
    figure('Name','P10 – Disturbances (Offset-free)');
    subplot(2,1,1); hold on; grid on;
    plot(t, d_true(1,:), 'k--','LineWidth',1.4,'DisplayName','True \DeltaF_3');
    plot(t, d_hat(1,:),  'b',  'LineWidth',1.5,'DisplayName','\hat{\DeltaF}_3');
    xline(t(k_step),'r--','d-step');
    ylabel('\DeltaF_3 [cm^3/s]');
    legend('Location','best');
    title(sprintf('F_3 disturbance step = %.1f cm^3/s', d_step_val(1)));

    subplot(2,1,2); hold on; grid on;
    plot(t, d_true(2,:), 'k--','LineWidth',1.4,'DisplayName','True \DeltaF_4');
    plot(t, d_hat(2,:),  'b',  'LineWidth',1.5,'DisplayName','\hat{\DeltaF}_4');
    xline(t(k_step),'r--','d-step');
    ylabel('\DeltaF_4 [cm^3/s]'); xlabel('Time [s]');
    legend('Location','best');
    title('F_4 disturbance (here no step)');
end

% -------------------------------------------------------------------------
function [h1_ref, h2_ref] = generate_stair_references( ...
                                    t, t_change_h1, values_h1, ...
                                       t_change_h2, values_h2)
% Simple staircase generator (same style as your previous problems)

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
