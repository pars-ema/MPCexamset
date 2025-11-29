function [x, y_true_abs, y_meas_abs, y_meas_dev] = simulate_sde_2( ...
        t, x0, u, d, p, sigma_meas, output_index, xs)
% SIMULATE_SDE_2
% Stochastic nonlinear Four-Tank System (Euler–Maruyama) with
% deviation measurements for Kalman filtering.
%
% Inputs:
%   t           : time vector [1×N]
%   x0          : initial state [4×1] (masses)
%   u           : inputs [2×1] or [2×N]     (F1,F2)
%   d           : disturbances [2×1] or [2×N] (F3,F4 mean)
%   p           : parameter vector/struct
%   sigma_meas  : measurement noise std for EACH TANK [4×1]
%   output_index: indices of measured heights (e.g. [1 2])
%   xs          : steady-state state [4×1] (for deviation)
%
% Outputs:
%   x           : true states [4×N] (masses)
%   y_true_abs  : true heights [4×N]
%   y_meas_abs  : measured heights (subset) [ny_meas×N]
%   y_meas_dev  : deviation heights (subset) [ny_meas×N],
%                 i.e. y_meas_abs - y_s_meas

    % ----------------------------------------------------
    % Basic setup
    % ----------------------------------------------------
    Ts = t(2) - t(1);
    N  = length(t);
    nx = 4;
    ny = 4;
    ny_meas = length(output_index);

    % Ensure u and d are 2×N
    if size(u,2) == 1
        u = repmat(u,1,N);
    end
    if size(d,2) == 1
        d = repmat(d,1,N);
    end

    % Preallocate
    x          = zeros(nx,N);
    y_true_abs = zeros(ny,N);
    y_meas_abs = zeros(ny_meas,N);
    y_meas_dev = zeros(ny_meas,N);

    % Initial state
    x(:,1) = x0;

    % ----------------------------------------------------
    % Disturbance Wiener increments for F3, F4
    % ----------------------------------------------------
    tf = t(end);
    [~,~,dW3] = ScalarStdWienerProcess(tf, N, 1);
    [~,~,dW4] = ScalarStdWienerProcess(tf, N, 1);
    sigma3 = 2;
    sigma4 = 2;

    % ----------------------------------------------------
    % Steady-state output and its measured part
    % ----------------------------------------------------
    y_s_abs  = FourTankSystemSensor(xs,p);   % 4×1 heights at steady state
    y_s_meas = y_s_abs(output_index);        % subset for measured heights

    % ----------------------------------------------------
    % Initial outputs
    % ----------------------------------------------------
    y_true_abs(:,1) = FourTankSystemSensor(x(:,1),p);

    y_meas_abs(:,1) = y_true_abs(output_index,1) + ...
        sigma_meas(output_index).*randn(ny_meas,1);

    % deviation = measured - steady-state measured
    y_meas_dev(:,1) = y_meas_abs(:,1) - y_s_meas;

    % ----------------------------------------------------
    % Euler–Maruyama integration loop
    % ----------------------------------------------------
    for k = 1:N-1

        % Nonlinear drift (deterministic part)
        f = FourTankSystemModified(0, x(:,k), u(:,k), d(:,k), p);

        % Process noise only on tank 3 and 4 (example choice)
        noise = [0;
                 0;
                 sigma3*dW3(k);
                 sigma4*dW4(k)];

        % State update
        x(:,k+1) = x(:,k) + Ts*f + noise;

        % True outputs (heights, all 4 tanks)
        y_true_abs(:,k+1) = FourTankSystemSensor(x(:,k+1),p);

        % Measured outputs (subset of heights)
        y_meas_abs(:,k+1) = y_true_abs(output_index,k+1) + ...
            sigma_meas(output_index).*randn(ny_meas,1);

        % Deviation measurements for KF
        y_meas_dev(:,k+1) = y_meas_abs(:,k+1) - y_s_meas;
    end

end
