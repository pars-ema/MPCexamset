clear; clc; close all;

% --- Load discrete-time models from your previous work ---
load('discrete_model_from_step.mat',          'Ad','Bd','Cd','Dd');  Ad_step=Ad; Bd_step=Bd; Cd_step=Cd; Dd_step=Dd;
load('discrete_model_from_linearization.mat', 'Ad','Bd','Cd','Dd');  Ad_lin =Ad; Bd_lin =Bd; Cd_lin =Cd; Dd_lin =Dd;

% Ts = 10;              
% nxS = size(Ad_step,1); nyS = size(Cd_step,1); nuS = size(Bd_step,2);
% nxL = size(Ad_lin ,1); nyL = size(Cd_lin ,1); nuL = size(Bd_lin ,2);
% N   = 800;                     % simulation length


%%
%% -------------------------------------------------------------
%  Setup: dimensions, noise covariances, simulation length
% --------------------------------------------------------------
Ts = 10;      % sampling time [s]
N  = 800;     % number of simulation steps

% Choose the model you want to simulate:
Ad = Ad_step;   Bd = Bd_step;   Cd = Cd_step;   Dd = Dd_step;
%Ad = Ad_lin;  Bd = Bd_lin;    Cd = Cd_lin;    Dd = Dd_lin;

nx = size(Ad,1);
ny = size(Cd,1);
nu = size(Bd,2);

% Disturbance/measurement noise covariance (tune these!)
Q = 1e-4 * eye(nx);     % process noise covariance
R = 1e-3 * eye(ny);     % measurement noise covariance
S = zeros(nx,ny);       % cross covariance (usually zero)

% Storage
x      = zeros(nx,N);     % true state
x_non_linear = zeros(nx,N);
y      = zeros(ny,N);     % measurements
u      = ones(nu,N).*250;     % inputs
xhat   = zeros(nx,N);     % KF estimate
Pstore = zeros(nx,nx,N);  % KF covariance history

% Initial values
x(:,1)    = zeros(nx,1);
xhat(:,1) = zeros(nx,1);
P         = 1*eye(nx);


%% -------------------------------------------------------------
%   Dynamic Kalman Filter loop
% --------------------------------------------------------------
for k = 1:N-1

    %% ---- Simulate true system (with stochastic disturbances) ----
    wk = mvnrnd(zeros(nx,1), Q)';   % process noise
    vk = mvnrnd(zeros(ny,1), R)';   % measurement noise

    y(:,k) = Cd * x(:,k) + vk;
    x(:,k+1) = Ad*x(:,k) + Bd*u(:,k) + wk;


    %% ---- Prediction (time update) ----
    xhat_pred = Ad*xhat(:,k) + Bd*u(:,k);
    P_pred    = Ad*P*Ad' + Q;   % if S = 0


    %% ---- Measurement Update ----
    ek   = y(:,k) - Cd*xhat_pred;                    % innovation
    Re   = Cd*P_pred*Cd' + R;                       % innovation covariance
    Kx   = P_pred*Cd'/Re;                           % state gain (Kf_x)
                                                   % (scalar inversion is ok)
    % Kf_w not needed unless you estimate disturbances


    xhat(:,k) = xhat_pred + Kx*ek;                  % updated estimate
    P         = P_pred - Kx*Re*Kx';                 % cov update

    Pstore(:,:,k) = P;                              % store
end

y(:,N)    = Cd*x(:,N);   % last measurement
xhat(:,N) = xhat(:,N-1); % extend


%% -------------------------------------------------------------
%  Static Kalman Filter gain (solve DARE)
% --------------------------------------------------------------
[Pinf,~,~] = dare(Ad', Cd', Q, R);
K_static = Pinf * Cd' / (Cd*Pinf*Cd' + R);

xhat_s = zeros(nx,N);
xhat_s(:,1) = zeros(nx,1);

for k = 1:N-1
    y(:,k) = Cd*x(:,k);   % reuse simulated measurement

    % Prediction
    xhat_pred = Ad*xhat_s(:,k) + Bd*u(:,k);

    % Update using constant gain
    ek_s = y(:,k) - Cd*xhat_pred;
    xhat_s(:,k+1) = xhat_pred + K_static * ek_s;
end


figure;
for i = 1:nx
    subplot(nx,1,i)
    plot(x(i,:), 'k', 'LineWidth',1.2); hold on;
    plot(xhat(i,:), 'r--', 'LineWidth',1.2);
    title(['State ' num2str(i)]);
    legend('True','KF Estimate');
end

%%
% -----------------------------------------------------------
% Parameters
% -----------------------------------------------------------
a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;        % [cm2]
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327; % [cm2]
g = 981; rho = 1.00;
gamma1 = 0.58; gamma2 = 0.72;
p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2; rho];
compare_linear_kf_with_sde(p)


function compare_linear_kf_with_sde(p)

    % -----------------------------------------------------------
    % 1. LOAD DISCRETE-TIME LINEAR MODEL (Problem 5)
    % -----------------------------------------------------------
    % We assume this model is: x_{k+1} = A x_k + B u_k,  y_k = C x_k
    % and that states are heights (or at least consistent with SDE model).
    %
    % From the slides, this corresponds to:
    %   x_{k+1} = A x_k + G w_k
    %   y_k     = C x_k + v_k
    % where we will add G and noise later.
    % -----------------------------------------------------------
    load('discrete_model_from_linearization.mat','Ad','Bd','Cd');

    Ad_lin = Ad;
    Bd_lin = Bd;
    
    % We only use measurements from tank 1 and 2:
    % C_meas x = [h1; h2]
    C_meas = Cd(1:2,:);   % 2 x nx
    ny = size(C_meas,1);

    % -----------------------------------------------------------
    % 2. SIMULATION PARAMETERS (NONLINEAR "TRUE" MODEL)
    % -----------------------------------------------------------
    t  = 0:10:1800;         % [s]
    Ts = 10;
    x0 = [100; 100; 100; 100];
    u  = [20; 20] * ones(1,length(t));  % constant pumps
    d  = [300; 300];                    % mean disturbances F3, F4

    sigma_meas = [5;5;5;5];             % measurement noise std for all 4 sensors

    % Nonlinear SDE simulation (truth)
    % simulate_sde implements the stochastic nonlinear model using Euler–Maruyama.
    [~, y_true, y_meas] = simulate_sde(t, x0, u, d, p, sigma_meas, false);

    % Extract the two measurements that the KF will actually use:
    y_meas_2 = y_meas(1:2,:);           % measured (noisy) h1, h2
    y_true_2 = y_true(1:2,:);           % true h1, h2

    % -----------------------------------------------------------
    % 3. DIMENSIONS
    % -----------------------------------------------------------
    nx = size(Ad_lin,1);

    % -----------------------------------------------------------
    % 4. NOISE MODELS (Q, R) AND G MATRIX
    % -----------------------------------------------------------
    % Measurement noise covariance R:
    % Slides: v_k ~ N(0, R)
    %
    % We only use h1, h2 -> use their variances:
    R = diag(sigma_meas(1:2).^2);   % 2 x 2

    % Process noise:
    % In the slides, x_{k+1} = A x_k + G w_k,  w_k ~ N(0, Q)
    % We do NOT try to perfectly reproduce the SDE noise here.
    % Instead, we treat process noise as generic modeling error + disturbances.
    %
    % Choice: G = I, Q = q * I
    % This matches the simple case in the slides (S = 0, G = I).
    % You can tune q to trade off smoothness vs. tracking.
    q = 1e-2;                    % <-- try 1e-3, 1e-2, 1e-1 to see effect
    G = eye(nx);
    Q = q * eye(nx);

    % -----------------------------------------------------------
    % 5. STATIC (TIME-INVARIANT) KALMAN FILTER (Slides 28–31)
    % -----------------------------------------------------------
    % Slides: Discrete Algebraic Riccati Equation (DARE)
    %
    %   P = A P A' + G Q G' - (A P C' + G S)(C P C' + R)^{-1}(A P C' + G S)'
    %
    % With S = 0, G = I, this simplifies and MATLAB dare uses:
    %   P = dare(A', C', G Q G', R)
    %
    % Then the stationary Kalman gain is:
    %   Re = C P C' + R
    %   K  = P C' Re^{-1}
    % -----------------------------------------------------------
    S = zeros(nx,ny);  %#ok<NASGU>  % no cross-covariance
    [P_ss,~,~] = dare(Ad_lin', C_meas', G*Q*G', R);   % P_ss = steady-state P

    Re_ss = C_meas * P_ss * C_meas' + R;             % Re = CPC' + R
    K_ss  = P_ss * C_meas' / Re_ss;                  % Kfx = PC' Re^{-1}

    % -----------------------------------------------------------
    % 6. DYNAMIC (TIME-VARYING) KALMAN FILTER (Slides 19–23)
    % -----------------------------------------------------------
    % Time-varying KF starts with some P0 and recursively updates:
    %
    % Measurement update:
    %   Re,k = C P_{k|k-1} C' + R
    %   K_k  = P_{k|k-1} C' Re,k^{-1}
    %   e_k  = y_k - C x̂_{k|k-1}
    %   x̂_{k|k} = x̂_{k|k-1} + K_k e_k
    %   P_{k|k} = P_{k|k-1} - K_k Re,k K_k'
    %
    % Time update:
    %   x̂_{k+1|k} = A x̂_{k|k} + B u_k
    %   P_{k+1|k} = A P_{k|k} A' + G Q G'
    % -----------------------------------------------------------
    xhat_dyn = zeros(nx,length(t));   % dynamic KF estimate
    P_dyn    = 10*eye(nx);            % initial covariance (tunable)

    xhat_dyn(:,1) = x0;

    % -----------------------------------------------------------
    % 7. STATIC KF ESTIMATE STORAGE
    % -----------------------------------------------------------
    xhat_stat = zeros(nx,length(t));
    xhat_stat(:,1) = x0;

    % -----------------------------------------------------------
    % 8. MAIN FILTERING LOOP
    % -----------------------------------------------------------
    for k = 1:length(t)-1

        % ------------------------- STATIC KF (steady-state K_ss) -------------------------
        % Time update:
        %   x̂_{k+1|k} = A x̂_{k|k} + B u_k
        x_pred_stat = Ad_lin * xhat_stat(:,k) + Bd_lin * u(:,k);

        % Measurement update with constant K_ss:
        %   e_k = y_k - C x̂_{k|k-1}
        yk = y_meas_2(:,k);
        e_stat = yk - C_meas * x_pred_stat;

        %   x̂_{k+1|k+1} = x̂_{k+1|k} + K_ss e_k
        xhat_stat(:,k+1) = x_pred_stat + K_ss * e_stat;

        % ------------------------- DYNAMIC KF (time-varying K_k) ------------------------
        % Time update:
        x_pred_dyn = Ad_lin * xhat_dyn(:,k) + Bd_lin * u(:,k);
        P_pred_dyn = Ad_lin * P_dyn * Ad_lin' + G * Q * G';

        % Measurement update as in slides 19–23:
        yk = y_meas_2(:,k);
        e_dyn = yk - C_meas * x_pred_dyn;               % ek = yk - Cx̂k|k-1

        Re_k = C_meas * P_pred_dyn * C_meas' + R;       % Re,k = C P C' + R
        K_k  = P_pred_dyn * C_meas' / Re_k;             % Kfx,k = P C' Re^{-1}

        xhat_dyn(:,k+1) = x_pred_dyn + K_k * e_dyn;     % x̂k|k = x̂k|k-1 + K ek
        P_dyn = P_pred_dyn - K_k * Re_k * K_k';         % Pk|k = Pk|k-1 - K Re K'
    end

    % -----------------------------------------------------------
    % 9. OUTPUTS FROM ESTIMATES
    % -----------------------------------------------------------
    y_est_stat = C_meas * xhat_stat;   % static KF output estimate (h1,h2)
    y_est_dyn  = C_meas * xhat_dyn;    % dynamic KF output estimate (h1,h2)

    % -----------------------------------------------------------
    % 10. PLOTS AND COMPARISONS
    % -----------------------------------------------------------
    figure;
    for i = 1:2
        subplot(2,2,i)
        hold on; grid on;
        plot(t, y_true_2(i,:), 'Color',[0.2 0.2 0.2], 'LineWidth', 1.5);  % true
        plot(t, y_meas_2(i,:), 'k:', 'LineWidth', 1);                    % noisy
        plot(t, y_est_stat(i,:), 'r--', 'LineWidth', 1.4);               % static KF
        plot(t, y_est_dyn(i,:),  'b-.', 'LineWidth', 1.4);               % dynamic KF
        xlabel('Time [s]');
        ylabel(sprintf('h_%d [cm]', i));
        title(sprintf('Tank %d: True vs Measured vs KF', i));
        legend('True (nonlinear SDE)','Measured (noisy)','Static KF','Dynamic KF', ...
               'Location','best');
    end

    % Innovation plots (to see if filter is "consistent")
    for i = 1:2
        subplot(2,2,2+i)
        hold on; grid on;
        innov_stat = y_meas_2(i,:) - y_est_stat(i,:);
        innov_dyn  = y_meas_2(i,:) - y_est_dyn(i,:);
        plot(t, innov_stat, 'r--', 'LineWidth', 1);
        plot(t, innov_dyn,  'b-.', 'LineWidth', 1);
        xlabel('Time [s]');
        ylabel(sprintf('Innovation h_%d', i));
        title(sprintf('Innovation for Tank %d', i));
        legend('Static KF innovation','Dynamic KF innovation','Location','best');
    end

end
