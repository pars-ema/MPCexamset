%% ========================================================================
% P12.3 – Prediction-Error Identification (PEM) for modified four-tank system
% ========================================================================

clear; clc; close all;

%% ----------------------------- Settings ---------------------------------
Ts   = 10;       % sample time [s]
Ns   = 250;      % number of samples
dt   = 0.05;     % internal integration step [s] for continuous propagation

nx = 4; ny = 2; nu = 2; nd = 2;

rng(1);

%% ------------------------ Physical parameters ---------------------------
% Same parameter vector layout as P8:
% p = [a1;a2;a3;a4;A1;A2;A3;A4;g;gamma1;gamma2;rho];

a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;        % [cm^2]
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327;% [cm^2]
g  = 981;                                                   % [cm/s^2]
rho = 1;                                                    % [g/cm^3]
gamma1 = 0.58;
gamma2 = 0.72;

p_true = [a1;a2;a3;a4;A1;A2;A3;A4;g;gamma1;gamma2;rho];

% Operating point (same style as P8)
u_s = [300; 300];     % [cm^3/s] (F1,F2)
d_s = [250; 250];     % [cm^3/s] (F3,F4) for nonlinear model

% Find steady-state masses for initial condition
xs_guess = 500*ones(nx,1);
opts = optimoptions('fsolve','Display','none');
x_s = fsolve(@(x)FourTankSystemModified(0,x,u_s,d_s,p_true), xs_guess, opts);

fprintf("Steady-state x_s (masses):\n"); disp(x_s);

% Measurement is lower tank heights (1 and 2)
h_s = FourTankSystemSensor(x_s, p_true);   % 4x1
y_s = h_s([1 2]);

%% ------------------------ Noise (stochastic sim) ------------------------
% Measurement noise covariance
R = diag([4 4]);    % heights noise

% Process noise: we inject diffusion on the disturbance inflow channels
% so that masses 3 and 4 get random inflow-like perturbations.
% sigma(x,u,p) dW  ~  G * L * dW, with G = [0 0; 0 0; 1 0; 0 1]
G = [0 0; 0 0; 1 0; 0 1];   % affects m3 and m4

% Choose diffusion intensity on the two channels
Qc = diag([25 25]);   % continuous-time intensity
Lw = chol(Qc, 'lower');

%% ------------------------ Input excitation ------------------------------
U = generate_prbs_inputs(Ns, nu, 'umin', [200;200], 'umax', [400;400], 'hold', 8);
Ddev = zeros(nd, Ns);  % keep disturbances at operating point in mean

%% ------------------------ Generate dataset ------------------------------
x0 = x_s + [50; -30; 20; -10];  % perturbed start for richer transient data

data = simulate_true_sde_dataset(x0, U, Ddev, Ts, dt, p_true, u_s, d_s, R, G, Lw);

% data fields:
%  data.t  : 1xNs
%  data.u  : 2xNs (absolute)
%  data.y_meas  : 2xNs (noisy measurements)
%  data.y_true  : 2xNs (true outputs without noise)
%  data.x0 : 4x1

%% ------------------------ plot simulated data ---------------------------
% plot containing true states and measurements (with and without noise) and the inputs
figure;
subplot(3,1,1);
plot(data.t, data.y_meas(1,:), 'o', 'DisplayName','y_1 measurement (noisy)');
hold on;
plot(data.t, data.y_true(1,:), '-', 'DisplayName','y_1 true (no noise)');
ylabel('h_1 [cm]'); legend('Location','best');

subplot(3,1,2);
plot(data.t, data.y_meas(2,:), 'o', 'DisplayName','y_2 measurement (noisy)');
hold on;
plot(data.t, data.y_true(2,:), '-', 'DisplayName','y_2 true (no noise)');
ylabel('h_2 [cm]'); legend('Location','best');  

subplot(3,1,3);
plot(data.t, data.u(1,:), '-', 'DisplayName','u_1 (F1)');
hold on;
plot(data.t, data.u(2,:), '-', 'DisplayName','u_2 (F2)');
ylabel('Inflow [cm^3/s]'); xlabel('time [s]'); legend('Location','best');

sgtitle('Simulated data for PEM identification');


%% ------------------------ Apply EKF for state estimation ------------------------
%apply EKF and plot results
p_theta = p_true; %use true parameters for EKF
ekf_data = run_ekf_on_simdata(data, p_theta, Ts, dt, R, Qc, G, d_s);

%plot EKF results
figure;
for i =1:2
    subplot(2,1,i); hold all;
    plot(data.t, data.y_meas(i,:), 'o', 'DisplayName','y_1 measurement (noisy)');
    plot(data.t, data.y_true(i,:), 'k-', 'DisplayName','y_1 true (no noise)');
    plot(data.t, ekf_data.y_hat(i,:), '-', 'DisplayName','y_1 EKF filter pred (true p)');
    plot(data.t, ekf_data.y_pred(i,:), '--', 'DisplayName','y_1 EKF pred (true p)');
    ylabel(sprintf('h_%d [cm]', i));
    legend('Location','best');
end

%do simila experiment but with wrong parameters
p_theta_wrong = p_true;
p_theta_wrong(1:4) = p_true(1:4) * 1.5; %use wrong parameters for EKF
ekf_data_wrong = run_ekf_on_simdata(data, p_theta_wrong, Ts, dt, R, Qc, G, d_s);
%plot EKF results with wrong parameters
figure;
for i =1:2
    subplot(2,1,i); hold all;
    plot(data.t, data.y_meas(i,:), 'o', 'DisplayName','y_1 measurement (noisy)');
    plot(data.t, data.y_true(i,:), 'k-', 'DisplayName','y_1 true (no noise)');
    plot(data.t, ekf_data_wrong.y_hat(i,:), '-', 'DisplayName','y_1 EKF filter pred (wrong p)');
    plot(data.t, ekf_data_wrong.y_pred(i,:), '--', 'DisplayName','y_1 EKF pred (wrong p)');
    ylabel(sprintf('h_%d [cm]', i));
    legend('Location','best');
end

%% ------------------------ Identify parameters via PEM ------------------------
% --- estimate ALL parameters ---
idx_est_all = 1:12;

p0=p_true;
% replace with perturbed initial guess
% 30% random perturbation
p0(idx_est_all) = p_true(idx_est_all) .* (1 + 0.30*randn(size(p_true(idx_est_all))));   
% avoid zeros/negatives
p0(idx_est_all) = max(p0(idx_est_all), 1e-6);           


% bounds must exist and be sensible for all parameters
% (example: +/- 5x around nominal)
lb = 0.2*p_true;  ub = 5.0*p_true;
% for gammas, force (0,1)
lb(10:11) = 0.01; ub(10:11) = 0.99;

[p_hat_all, pem_out_all] = pem_identify_params(data, p0, idx_est_all, lb, ub, Ts, dt, R, Qc, G, d_s);

disp(table(p_true(idx_est_all), p_hat_all(idx_est_all), ...
    (p_hat_all(idx_est_all)-p_true(idx_est_all))./p_true(idx_est_all), ...
    'VariableNames', {'p_true','p_hat','rel_err'}));

    ekf_ident_all = run_ekf_on_simdata(data, p_hat_all, Ts, dt, R, Qc, G, d_s);

figure;
for i = 1:2
    subplot(2,1,i); hold on; grid on;
    plot(data.t, data.y_meas(i,:), 'o', 'DisplayName',sprintf('y_%d meas (noisy)',i));
    plot(data.t, data.y_true(i,:), 'k-', 'DisplayName',sprintf('y_%d true',i));
    plot(data.t, ekf_ident_all.y_hat(i,:), '-',  'DisplayName',sprintf('y_%d EKF (filtered)',i));
    plot(data.t, ekf_ident_all.y_pred(i,:), '--','DisplayName',sprintf('y_%d EKF (one-step)',i));
    ylabel(sprintf('h_%d [cm]', i));
    legend('Location','best');
end
xlabel('time [s]');


figure('Name','PEM relative parameter updates');
theta_hist = pem_out_all.history; %nEst x nIters
true_theta = p_true(idx_est_all);
iters = 0:(size(theta_hist,2)-1);
param_names = {'a_1','a_2','a_3','a_4','A_1','A_2','A_3','A_4','g','\gamma_1','\gamma_2','rho'};
for j = idx_est_all
    subplot(4,3,j); hold on; grid on;
    est = theta_hist(j,:);
    plot(iters, est, '-');
    yline(true_theta(j), 'k--', 'LineWidth', 1.0);
    title(param_names{idx_est_all(j)});
    xlabel('iteration');
    ylabel('estimate');
    %set ylim to lb and ub relative
    ylim([lb(j), ub(j)]);
end
sgtitle('Parameter updates');

%%
nRuns = 1;
idx_est_all = 10:11;
nEst = numel(idx_est_all);

param_names = {'a_1','a_2','a_3','a_4','A_1','A_2','A_3','A_4','g','\gamma_1','\gamma_2','rho'};
true_theta = p_true(idx_est_all);

% bounds (same as you already do)
lb = 0.2*p_true;  ub = 5.0*p_true;
lb(10:11) = 0.01; ub(10:11) = 0.99;

% store histories (each run can have different iteration count)
theta_hist_runs = cell(nRuns,1);
p_hat_runs      = zeros(numel(p_true), nRuns);
resnorm_runs    = nan(nRuns,1);

for r = 1:nRuns
    rng(r); % different random initial point

    % random initial guess for the estimated subset
    p0 = p_true;
    p0(idx_est_all) = p_true(idx_est_all) .* (1 + 0.30*randn(nEst,1));
    p0(idx_est_all) = max(p0(idx_est_all), 1e-6);

    [p_hat, pem_out] = pem_identify_params(data, p0, idx_est_all, lb, ub, Ts, dt, R, Qc, G, d_s);

    p_hat_runs(:,r) = p_hat;
    theta_hist_runs{r} = pem_out.history;   % nEst x nIters
    if isfield(pem_out,'resnorm')
        resnorm_runs(r) = pem_out.resnorm;
    end

    fprintf('Run %d/%d done. Final est = %s\n', r, nRuns, mat2str(p_hat(idx_est_all)',4));
end

% -------- Plot: overlay parameter updates for all runs ----------
figure('Name','PEM parameter updates (10 runs)');
for j = 1:length(idx_est_all)
    subplot(4,3,j); hold on; grid on;

    % plot each run trajectory
    for r = 1:nRuns
        th = theta_hist_runs{r};              % nEst x nIters
        it = 0:(size(th,2)-1);
        plot(it, th(j,:), '-');              % overlay
    end

    % true parameter line
    yline(true_theta(j), 'k--', 'LineWidth', 1.0);

    title(param_names{idx_est_all(j)});
    xlabel('iteration'); ylabel('estimate');

    % use bounds for THIS parameter index (careful: lb/ub are full p)
    ylim([lb(idx_est_all(j)), ub(idx_est_all(j))]);
end
sgtitle(sprintf('Parameter updates over %d PEM runs (different initial guesses)', nRuns));




%% ========================================================================
% =============================== FUNCTIONS ================================
% ========================================================================

function data = simulate_true_sde_dataset(x0, U, Ddev, Ts, dt, p_true, u_s, d_s, R, G, Lw)
    % Euler–Maruyama simulation with ZOH inputs.
    % xdot = f(x,u,d,p) + G*Lw * w(t),  w(t) is white noise -> dW ~ sqrt(dt) N(0,I)

    Ns = size(U,2);
    x  = x0;
    y_true  = zeros(2,Ns);
    y_meas = zeros(2,Ns);
    t  = (0:Ns-1)*Ts;

    nSteps = ceil(Ts/dt);
    h = Ts/nSteps;

    for k = 1:Ns
        u_abs = U(:,k);
        d_abs = d_s + Ddev(:,k);

        for i = 1:nSteps
            dW = sqrt(h) * randn(2,1);              % 2-dim Brownian increment
            drift = FourTankSystemModified(0, x, u_abs, d_abs, p_true);
            x = x + h*drift + (G*Lw)*dW;
        end

        vk = chol(R,'lower') * randn(2,1);
        hk = FourTankSystemSensor(x, p_true);
        % y(:,k) = hk([1 2]) + vk;
        y_true(:,k) = hk([1 2]);
        y_meas(:,k) = hk([1 2]) + vk;
    end

    data.t  = t;
    data.u  = U;
    data.y_meas  = y_meas;
    data.y_true  = y_true;
    data.x0 = x0;
end

function U = generate_prbs_inputs(Ns, nu, varargin)
    % Simple random-step PRBS-like inputs with hold time
    p = inputParser;
    addParameter(p,'umin', zeros(nu,1));
    addParameter(p,'umax', ones(nu,1));
    addParameter(p,'hold', 10);
    parse(p,varargin{:});
    umin = p.Results.umin(:); umax = p.Results.umax(:);
    holdN = p.Results.hold;

    U = zeros(nu,Ns);
    k = 1;
    while k <= Ns
        u = umin + (umax-umin).*rand(nu,1);
        k2 = min(Ns, k+holdN-1);
        U(:,k:k2) = repmat(u, 1, k2-k+1);
        k = k2+1;
    end
end

% Extended kalman filter (EKF) for state estimation
function ekf = run_ekf_on_simdata(data, p_hat, Ts, dt, R, Qc, G, d_s)
    % Continuous–discrete EKF on simulated dataset.
    % Returns both:
    %   y_pred(:,k) = h(x_{k|k-1})  (one-step-ahead prediction)
    %   y_hat(:,k)  = h(x_{k|k})    (filtered estimate after update)

    Ns = numel(data.t);
    nx = 4; ny = 2;

    % Constant measurement Jacobian for y=[h1;h2] with h_i = m_i/(rho*A_i)
    rho = p_hat(12);
    A1  = p_hat(5);
    A2  = p_hat(6);
    H = [1/(rho*A1), 0, 0, 0;
        0, 1/(rho*A2), 0, 0];

    % Continuous-time process noise covariance in state coordinates
    Qekf = G * Qc * G';   % (4x4)

    % Preallocate outputs
    xhat  = zeros(nx,Ns);
    xpred = zeros(nx,Ns);
    P     = zeros(nx,nx,Ns);

    y_pred = zeros(ny,Ns);   % one-step-ahead predictions
    y_hat  = zeros(ny,Ns);   % filtered output estimates

    innov = zeros(ny,Ns);
    Sall  = zeros(ny,ny,Ns);
    Kall  = zeros(nx,ny,Ns);

    % Init
    xk = data.x0(:);
    Pk = 1e2*eye(nx);

    % Prediction substeps
    nSteps = max(1, ceil(Ts/dt));
    h = Ts/nSteps;

    for k = 1:Ns
        uk = data.u(:,k);
        y_meas = data.y_meas(:,k);

        % ------------------- Prediction: propagate mean + covariance -----------
        x = xk;
        Pk_pred = Pk;

        for i = 1:nSteps
            % Mean propagation
            f = FourTankSystemModified(0, x, uk, d_s, p_hat);
            x = x + h*f;
            x = max(x, 1e-12);

            % Linearization A(t) = df/dx (numerical finite differences)
            %Acont = numerical_jacobian(@(xx) FourTankSystemModified(0, xx, uk, d_s, p_hat), x);
            Acont = FourTankSystemModified_Jacobian_f(x, p_hat);

            % Covariance propagation: Pdot = A P + P A' + Q
            Pk_pred = Pk_pred + h*(Acont*Pk_pred + Pk_pred*Acont' + Qekf);
            Pk_pred = 0.5*(Pk_pred + Pk_pred');  % symmetrize
        end

        % Predicted state
        xpred(:,k) = x;

        % One-step-ahead predicted measurement y_{k|k-1}
        hk = FourTankSystemSensor(x, p_hat);       % 4x1 heights
        yk_pred = hk([1 2]);                       % 2x1
        y_pred(:,k) = yk_pred;

        % Innovation (residual)
        ek = y_meas - yk_pred;
        innov(:,k) = ek;

        % ------------------- Update: discrete measurement ----------------------
        S = H*Pk_pred*H' + R;
        K = (Pk_pred*H')/S;

        % Updated state estimate
        xk = x + K*ek;

        % Joseph covariance update (robust)
        I = eye(nx);
        Pk = (I - K*H)*Pk_pred*(I - K*H)' + K*R*K';
        Pk = 0.5*(Pk + Pk');

        % Store filtered
        xhat(:,k) = xk;
        P(:,:,k)  = Pk;
        Sall(:,:,k) = S;
        Kall(:,:,k) = K;

        % Filtered output estimate y_{k|k}
        hk_f = FourTankSystemSensor(xk, p_hat);
        y_hat(:,k) = hk_f([1 2]);
    end

    % pack outputs
    ekf.xhat   = xhat;
    ekf.xpred  = xpred;
    ekf.P      = P;

    ekf.y_pred = y_pred;
    ekf.y_hat  = y_hat;

    ekf.innov  = innov;
    ekf.S      = Sall;
    ekf.K      = Kall;

end



function A = FourTankSystemModified_Jacobian_f(x,p)
    % Analytical Jacobian A = df/dx for modified four-tank mass model, x=[m1;m2;m3;m4]
    % p=[a1;a2;a3;a4; A1;A2;A3;A4; g; gamma1;gamma2; rho]

    m  = max(x(3:4),1e-12); % m3, m4 (avoid sqrt/div-by-zero)
    a  = p(3:4);  % a3, a4 
    At = p(7:8);  % A3, A4
    g = p(9); 
    rho = p(12);

    dq  = rho .* a .* sqrt(2*g./(rho.*At)) ./ (2*sqrt(m)); % dq_i/dm_i

    A = zeros(4); % df/dx
    A(1,3) = dq(1); % df1/dm3
    A(2,4) = dq(2); % df2/dm4
end

function H = FourTankSystemModified_Jacobian_h(x,p)
    % Analytical Jacobian H = dh/dx for modified four-tank mass model, x=[m1;m2;m3;m4]
    % p=[a1;a2;a3;a4; A1;A2;A3;A4; g; gamma1;gamma2; rho]

    At = p(5:6);  % A1, A2
    rho = p(12);

    dh = 1./(rho*At); % dh_i/dm_i

    A = zeros(2,4); % dh/dx
    A(1,1) = dh(1); % dh1/dm1
    A(1,2) = dh(2); % df2/dm2
end

function A = numerical_jacobian(fun, x)
    %NUMERICAL_JACOBIAN  Finite-difference Jacobian of fun(x) w.r.t. x.
    % fun: R^n -> R^m  (here m=n=4 for your dynamics)
    % x  : n-by-1
    
    fx = fun(x);
    n  = numel(x);
    m  = numel(fx);
    A  = zeros(m,n);
    
    eps0 = 1e-6;
    for i = 1:n
        dx = zeros(n,1);
        hi = eps0*(1 + abs(x(i)));   % scale step with magnitude
        dx(i) = hi;
        A(:,i) = (fun(x + dx) - fx)/hi;
    end
end


function [p_hat, out] = pem_identify_params(data, p0, idx_est, lb, ub, Ts, dt, R, Qc, G, d_s, opts)

    if nargin < 13 || isempty(opts)
        opts = optimoptions('lsqnonlin','Display','iter', ...
            'MaxFunctionEvaluations', 400, ...
            'FunctionTolerance', 1e-8, ...
            'StepTolerance', 1e-8);
    end

    theta0 = p0(idx_est);
    lb_th  = lb(idx_est);
    ub_th  = ub(idx_est);

    W = chol(inv(R), 'lower');  % OK (you can also do W = inv(chol(R,'lower')))

    resfun = @(theta) pem_residuals(theta, p0, idx_est, data, Ts, dt, R, Qc, G, d_s, W);

    % --------- capture iteration history via OutputFcn (version-robust) ----------
    theta_hist = [];

    function stop = outfun(x, optimValues, state) %#ok<INUSD>
        stop = false;

        if strcmp(state,'init')
            theta_hist = x(:);
        elseif strcmp(state,'iter')
            theta_hist(:,end+1) = x(:);
        end
    end

    opts = optimoptions(opts, 'OutputFcn', @outfun);

    % Run optimizer
    [theta_hat, resnorm, residual, exitflag, output] = lsqnonlin(resfun, theta0, lb_th, ub_th, opts);

    % Unpack full parameter vector
    p_hat = p0;
    p_hat(idx_est) = theta_hat;

    % Output
    out.resnorm  = resnorm;
    out.residual = residual;
    out.exitflag = exitflag;
    out.output   = output;     % contains iterations, funcCount, message, etc.
    out.history  = theta_hist;       % contains per-iteration updates
end


% ---------- helper: PEM residuals ----------
function r = pem_residuals(theta, p_base, idx_est, data, Ts, dt, R, Qc, G, d_s, W)
    p = p_base;
    p(idx_est) = theta;

    % One-step-ahead predictions from EKF
    ekf = run_ekf_on_simdata(data, p, Ts, dt, R, Qc, G, d_s);

    % Prediction errors (innovations) based on ONE-STEP predictor y_pred
    E = data.y_meas - ekf.y_pred;   % (2 x Ns)

    % Weight/whiten each error vector with R^{-1/2}
    Ew = W * E;

    % Stack into a long vector for lsqnonlin
    r = Ew(:);

end

% ---------- helper transforms for fallback ----------
function y = logit(z)
    y = 1./(1+exp(-z));
end

function z = invlogit(x)
    x = min(max(x,1e-9),1-1e-9);
    z = log(x./(1-x));
end
