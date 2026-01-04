%% ========================================================================
% P12 – Bound-constrained Nonlinear MPC (NMPC) for modified four-tank system
% ========================================================================

clear; clc; close all;

set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesColor','w');
set(groot,'defaultAxesXColor','k');
set(groot,'defaultAxesYColor','k');
set(groot,'defaultTextColor','k');
set(groot,'defaultLegendTextColor','k');
set(groot,'defaultLegendColor','w');

%% ----------------------------- Settings ---------------------------------
Ts   = 10;        % sample time [s]
Tsim = 300;       % number of MPC steps
t    = (0:Tsim-1)*Ts;

% Large dt makes it faster (coarser integration)
dt_plant = 1;     % integration step for plant simulation [s]
dt_ekf   = 1;     % integration step inside EKF prediction [s]
dt_nmpc  = 1;     % integration step inside NMPC prediction [s]

rng(0);           % reproducibility

nx = 4; ny = 2; nu = 2; nd = 2;

%% ------------------------ Physical parameters ---------------------------
a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;         % [cm^2]
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327; % [cm^2]
g  = 981;                                                    % [cm/s^2]
rho = 1;                                                     % [g/cm^3]
gamma1 = 0.58;
gamma2 = 0.72;

p = [a1;a2;a3;a4;A1;A2;A3;A4;g;gamma1;gamma2;rho];

% Operating point
u_s = [300;300];     % [cm^3/s] (F1, F2)
d_s = [250;250];     % [cm^3/s] (F3, F4)

% Find steady-state masses xs
xs_guess = 500*ones(nx,1);
opts = optimoptions('fsolve','Display','none');
xs = fsolve(@(x)FourTankSystemModified(0,x,u_s,d_s,p), xs_guess, opts);

fprintf("Steady-state xs (masses):\n"); disp(xs);

% Steady-state measured heights z_s = [h1; h2]
h_s = FourTankSystemSensor(xs,p);
zs  = h_s([1 2]);

fprintf("Steady-state zs = [h1_s; h2_s]:\n"); disp(zs);

%% ------------------------ NMPC tuning -----------------------------------
N = 15;                     % horizon

Q   = diag([1 1]);           % output tracking weight (h1,h2)
R   = 1e-4*eye(nu);          % input magnitude (u - us)
Rdu = 1e-2*eye(nu);          % input move penalty (delta u)

u_min_abs = [0;0];
u_max_abs = [600;600];

opt_nlp = optimoptions('fmincon', ...
    'Display','none', ...
    'Algorithm','sqp', ...
    'MaxIterations', 30, ...
    'MaxFunctionEvaluations', 8000, ...
    'StepTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-6);

%% ------------------------ EKF tuning (offset-free) ----------------------
R_meas = diag([4 4]);        % measurement noise covariance for [h1; h2]

Qx_ekf = 1e-2*eye(nx);       % process noise on x  (tune)
Qd_ekf = 0.5*eye(nd);        % random-walk noise on d (TUNE: larger => faster offset adaptation)

%% ========================================================================
% Experiment A: Reference steps, disturbances = 0
% ========================================================================

disp('=== NMPC – EXPERIMENT A: Reference steps, disturbances = 0 ===');

t_change_h1 = [0 600 1600];
values_h1   = [zs(1) zs(1)*1.3 zs(1)*0.8];

t_change_h2 = [0 600 1300];
values_h2   = [zs(2) zs(2)*1.3 zs(2)*1.0];

[h1_ref_A, h2_ref_A] = generate_stair_references( ...
    t, t_change_h1, values_h1, ...
       t_change_h2, values_h2);

z_ref_abs_A = [h1_ref_A; h2_ref_A];        % absolute ref
d_dev_A     = zeros(nd,Tsim);              % no disturbance dev

simA = simulateNMPC_nonlinear_offsetfree( ...
    xs, u_s, d_s, p, ...
    z_ref_abs_A, d_dev_A, ...
    Ts, dt_plant, dt_ekf, dt_nmpc, ...
    N, Q, R, Rdu, ...
    u_min_abs, u_max_abs, ...
    R_meas, Qx_ekf, Qd_ekf, opt_nlp);

%% ========================================================================
% Experiment B: Same references + disturbance step later
% ========================================================================

disp('=== NMPC – EXPERIMENT B: Reference steps + disturbance step ===');

z_ref_abs_B = z_ref_abs_A;

d_dev_B = zeros(nd,Tsim);
k_step = round(2000/Ts);                 % time index of disturbance step
d_step_val = [50; -30];                  % delta disturbance on [F3;F4]
d_dev_B(:,k_step:end) = repmat(d_step_val,1,Tsim-k_step+1);

simB = simulateNMPC_nonlinear_offsetfree( ...
    xs, u_s, d_s, p, ...
    z_ref_abs_B, d_dev_B, ...
    Ts, dt_plant, dt_ekf, dt_nmpc, ...
    N, Q, R, Rdu, ...
    u_min_abs, u_max_abs, ...
    R_meas, Qx_ekf, Qd_ekf, opt_nlp);

%% ========================================================================
% PLOTS
% ========================================================================

% ---------- Experiment A ----------
fig = figure;
t = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

figure('Name','NMPC_ExpA_outputs');
sgtitle('NMPC (offset-free) – Experiment A – Reference steps, d = 0');

subplot(3,1,1); hold on; grid on;
plot(t, simA.y_meas(1,:), 'kx','DisplayName','h_1 meas');
plot(t, z_ref_abs_A(1,:), 'k--','LineWidth',1.2, 'DisplayName','h_1 ref');
plot(t, simA.y_abs(1,:),  'b',  'LineWidth',1.4, 'DisplayName','h_1');
ylabel('h_1 [cm]');
legend('Location','best');

subplot(3,1,2); hold on; grid on;
plot(t, simA.y_meas(2,:), 'kx','DisplayName','h_2 meas');
plot(t, z_ref_abs_A(2,:), 'k--','LineWidth',1.2, 'DisplayName','h_2 ref');
plot(t, simA.y_abs(2,:),  'b',  'LineWidth',1.4, 'DisplayName','h_2');
ylabel('h_2 [cm]');
legend('Location','best');

subplot(3,1,3); hold on; grid on;
plot(t, simA.u(1,:), 'b','LineWidth',1.4,'DisplayName','u_1');
plot(t, simA.u(2,:), 'r','LineWidth',1.4,'DisplayName','u_2');
yline(u_s(1),'k--','LineWidth',1.0);
yline(u_s(2),'k--','LineWidth',1.0);
yline(0,'k:');
yline(600,'k:');
ylabel('u [cm^3/s]'); xlabel('Time [s]');
legend('Location','best');
% save in figures folder
exportgraphics(t, 'NMPC_ExpA_outputs.pdf', 'ContentType','vector');

t = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
figure('Name','NMPC_ExpA_dhat');
sgtitle('NMPC (offset-free) – Experiment A – Estimated disturbances');
subplot(2,1,1); hold on; grid on;
plot(t, simA.d_hat(1,:), 'LineWidth',1.4);
yline(d_s(1),'k--','LineWidth',1.0);
ylabel('\hat F_3 [cm^3/s]');
subplot(2,1,2); hold on; grid on;
plot(t, simA.d_hat(2,:), 'LineWidth',1.4);
yline(d_s(2),'k--','LineWidth',1.0);
ylabel('\hat F_4 [cm^3/s]'); xlabel('Time [s]');

% ---------- Experiment B ----------
figure('Name','NMPC_ExpB_outputs');
sgtitle('NMPC (offset-free) – Experiment B – Disturbance step');

subplot(3,1,1); hold on; grid on;
plot(t, simB.y_meas(1,:), 'kx','DisplayName','h_1 meas');
plot(t, z_ref_abs_B(1,:), 'k--','LineWidth',1.2,'DisplayName','h_1 ref');
plot(t, simB.y_abs(1,:),  'b',  'LineWidth',1.4,'DisplayName','h_1');
xline(t(k_step),'r--','Step in d','LabelVerticalAlignment','bottom');
ylabel('h_1 [cm]');
legend('Location','best');

subplot(3,1,2); hold on; grid on;
plot(t, simB.y_meas(2,:), 'kx','DisplayName','h_2 meas');
plot(t, z_ref_abs_B(2,:), 'k--','LineWidth',1.2,'DisplayName','h_2 ref');
plot(t, simB.y_abs(2,:),  'b',  'LineWidth',1.4,'DisplayName','h_2');
xline(t(k_step),'r--','Step in d','LabelVerticalAlignment','bottom');
ylabel('h_2 [cm]');
legend('Location','best');

subplot(3,1,3); hold on; grid on;
plot(t, simB.u(1,:), 'b','LineWidth',1.4,'DisplayName','u_1');
plot(t, simB.u(2,:), 'r','LineWidth',1.4,'DisplayName','u_2');
xline(t(k_step),'r--','Step in d');
yline(0,'k:'); yline(600,'k:');
ylabel('u [cm^3/s]'); xlabel('Time [s]');
legend('Location','best');
% save in figures folder
saveas(gcf, 'NMPC_ExpB_outputs.png');

figure('Name','NMPC_ExpB_dhat');
sgtitle('NMPC (offset-free) – Experiment B – True vs estimated disturbances');

subplot(2,1,1); hold on; grid on;
plot(t, d_s(1) + d_dev_B(1,:), 'k--','LineWidth',1.2,'DisplayName','True F_3');
plot(t, simB.d_hat(1,:),       'b',  'LineWidth',1.4,'DisplayName','\hat F_3');
xline(t(k_step),'r--','Step');
ylabel('F_3 [cm^3/s]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
plot(t, d_s(2) + d_dev_B(2,:), 'k--','LineWidth',1.2,'DisplayName','True F_4');
plot(t, simB.d_hat(2,:),       'b',  'LineWidth',1.4,'DisplayName','\hat F_4');
xline(t(k_step),'r--','Step');
ylabel('F_4 [cm^3/s]'); xlabel('Time [s]');
legend('Location','best');
saveas(gcf, 'P12_4_ExpB_dhat.png');

%% ========================================================================
% ============================= FUNCTIONS =================================
% ========================================================================

function sim = simulateNMPC_nonlinear_offsetfree( ...
    x0, u_s, d_s, p, ...
    z_ref_abs, d_dev, ...
    Ts, dt_plant, dt_ekf, dt_nmpc, ...
    N, Q, R, Rdu, ...
    umin, umax, ...
    R_meas, Qx_ekf, Qd_ekf, opt_nlp)

    nx = 4; ny = 2; nu = 2; nd = 2;
    nxe = nx + nd;
    Tsim = size(z_ref_abs,2);

    x_abs   = zeros(nx, Tsim+1);
    y_abs   = zeros(ny,Tsim);
    y_meas  = zeros(ny,Tsim);
    u_abs   = zeros(nu,Tsim);

    xehat   = zeros(nxe, Tsim+1);      % [x; d]
    Pe      = zeros(nxe,nxe,Tsim+1);
    d_hat   = zeros(nd,Tsim);

    % init
    x_abs(:,1)  = x0(:);

    xehat(:,1)  = [x0(:); d_s(:)];                 % start at nominal disturbance
    Pe(:,:,1)   = blkdiag(1e2*eye(nx), 1e2*eye(nd));

    Lr = chol(R_meas,'lower');

    % warm start
    Uwarm  = repmat(u_s(:), N, 1);
    u_prev = u_s(:);

    for k = 1:Tsim
        % ---- measurement ----
        hk = FourTankSystemSensor(x_abs(:,k), p);
        y_true = hk([1 2]);
        yk = y_true + Lr*randn(ny,1);

        y_abs(:,k)  = y_true;
        y_meas(:,k) = yk;

        % ---- EKF update (augmented, offset-free) ----
        [xe_k, Pe_k] = ekf_onestep_cd_aug(xehat(:,k), Pe(:,:,k), u_prev, yk, Ts, dt_ekf, p, R_meas, Qx_ekf, Qd_ekf);
        xehat(:,k) = xe_k;
        Pe(:,:,k)  = Pe_k;

        xhat_k = xe_k(1:nx);
        dhat_k = xe_k(nx+1:end);          % absolute [F3;F4] estimate
        d_hat(:,k) = dhat_k;

        % ---- reference over horizon ----
        zref = zeros(ny,N);
        for j = 1:N
            idx = min(k+j-1, Tsim);
            zref(:,j) = z_ref_abs(:,idx);
        end

        % ---- NMPC solve (absolute inputs) ----
        lb = repmat(umin(:), N, 1);
        ub = repmat(umax(:), N, 1);
        U0 = min(max(Uwarm, lb), ub);

        % IMPORTANT: predict with current disturbance estimate dhat_k
        Jfun = @(Uvec) nmpc_cost_nl(Uvec, xhat_k, u_prev, zref, Ts, dt_nmpc, dhat_k, p, Q, R, Rdu, u_s);

        Uopt = fmincon(Jfun, U0, [], [], [], [], lb, ub, [], opt_nlp);

        u_k = Uopt(1:nu);
        u_abs(:,k) = u_k;

        % warm start shift
        Uwarm  = [Uopt(nu+1:end); Uopt(end-nu+1:end)];
        u_prev = u_k;

        % ---- propagate TRUE nonlinear plant one sample ----
        d_k = d_s + d_dev(:,k);  % true plant disturbance (unknown to controller)
        x_abs(:,k+1) = plant_step_det(x_abs(:,k), u_k, d_k, Ts, dt_plant, p);

        % carry forward estimator states (next iter updates)
        xehat(:,k+1) = xehat(:,k);
        Pe(:,:,k+1)  = Pe(:,:,k);
    end

    sim = struct();
    sim.x_abs   = x_abs;
    sim.y_abs   = y_abs;
    sim.y_meas  = y_meas;
    sim.u       = u_abs;
    sim.d_hat   = d_hat;
end

function J = nmpc_cost_nl(Uvec, x0, u_prev, zref, Ts, dt, d_pred, p, Q, R, Rdu, u_s)
    nu = 2;
    N  = size(zref,2);
    U  = reshape(Uvec, nu, N);

    nSteps = max(1, ceil(Ts/dt));
    h = Ts/nSteps;

    x = x0(:);
    J = 0;

    for k = 1:N
        uk = U(:,k);

        for i = 1:nSteps
            f = FourTankSystemModified(0, x, uk, d_pred, p);
            x = max(x + h*f, 1e-12);
        end

        hk = FourTankSystemSensor(x, p);
        z  = hk([1 2]);

        e  = z - zref(:,k);
        ui = uk - u_s(:);

        if k == 1
            du = uk - u_prev(:);
        else
            du = uk - U(:,k-1);
        end

        J = J + e'*Q*e + ui'*R*ui + du'*Rdu*du;
    end
end

function x_next = plant_step_det(x, u, d, Ts, dt, p)
    nSteps = max(1, ceil(Ts/dt));
    h = Ts/nSteps;
    xk = x(:);
    for i = 1:nSteps
        f = FourTankSystemModified(0, xk, u, d, p);
        xk = max(xk + h*f, 1e-12);
    end
    x_next = xk;
end

function [xe_upd, Pe_upd] = ekf_onestep_cd_aug(xe_k, Pe_k, uk, yk, Ts, dt, p, R, Qx, Qd)
% Continuous–discrete EKF for augmented state xe=[x; d], where d=[F3;F4] is unknown.
% d is modeled as random walk (constant between samples with Qd diffusion).

    nx = 4; nd = 2; ny = 2;

    xk = xe_k(1:nx);
    dk = xe_k(nx+1:end);

    % Measurement Jacobian He = dh/d[x;d]
    rho = p(12); A1 = p(5); A2 = p(6);
    Hx = [1/(rho*A1), 0, 0, 0;
          0, 1/(rho*A2), 0, 0];
    He = [Hx, zeros(ny,nd)];

    Qe = blkdiag(Qx, Qd);

    nSteps = max(1, ceil(Ts/dt));
    h = Ts/nSteps;

    xpred = xk(:);
    dpred = dk(:);
    Ppred = Pe_k;

    for i = 1:nSteps
        % mean propagation (d constant between samples)
        f = FourTankSystemModified(0, xpred, uk, dpred, p);
        xpred = max(xpred + h*f, 1e-12);

        % Jacobians
        Axx = numerical_jacobian(@(xx) FourTankSystemModified(0, xx, uk, dpred, p), xpred);
        Axd = numerical_jacobian(@(dd) FourTankSystemModified(0, xpred, uk, dd, p), dpred);

        A = [Axx, Axd;
             zeros(nd,nx), zeros(nd,nd)];

        Ppred = Ppred + h*(A*Ppred + Ppred*A' + Qe);
        Ppred = 0.5*(Ppred + Ppred');
    end

    xepred = [xpred; dpred];

    % predicted measurement
    ypred = He*xepred;

    % update
    S = He*Ppred*He' + R;
    K = (Ppred*He')/S;

    xe_upd = xepred + K*(yk - ypred);

    I = eye(nx+nd);
    Pe_upd = (I - K*He)*Ppred*(I - K*He)' + K*R*K';
    Pe_upd = 0.5*(Pe_upd + Pe_upd');
end

function A = numerical_jacobian(fun, x)
    fx = fun(x);
    n  = numel(x);
    m  = numel(fx);
    A  = zeros(m,n);

    eps0 = 1e-6;
    for i = 1:n
        dx = zeros(n,1);
        hi = eps0*(1 + abs(x(i)));
        dx(i) = hi;
        A(:,i) = (fun(x + dx) - fx)/hi;
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
