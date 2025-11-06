function simulate_pid_stochastic(t, x0, u_s, p, h_ref_traj, Kp, Ki, Kd, outputFolder)
% SIMULATE_PID_STOCHASTIC
% PID control for nonlinear four-tank system with piecewise-constant stochastic disturbances.
%
% Inputs:
%   t            : time vector [s]
%   x0           : initial state [4x1]
%   u_s          : steady-state pump flows [2x1]
%   p            : parameter struct
%   h_ref_traj   : reference trajectory [2xN] (for tanks 1 & 2)
%   Kp, Ki, Kd   : PID gain vectors [2x1]
%   outputFolder : directory for saving figures
%
% The controller operates in deviation variables:
%   du = Kp.*e + Ki.*∫e dt + Kd.*de/dt
%   u  = u_s + du, with saturation limits [0, 1000].
%
% Author: MPC Project Team
% -------------------------------------------------------------

    if nargin < 9 || isempty(outputFolder)
        outputFolder = fullfile('figures','problem_3','pid_stochastic');
    end
    if ~exist(outputFolder,'dir'); mkdir(outputFolder); end

    % --- Time setup ---
    Ts = t(2) - t(1);
    N  = length(t);

    % --- Initialization ---
    nx = 4; ny = 4; nu = 2;
    x = zeros(nx,N);
    y_true = zeros(ny,N);
    y_meas = zeros(ny,N);
    u = zeros(nu,N);
    int_e = zeros(2,1);
    prev_e = zeros(2,1);
    x(:,1) = x0;

    % --- Measurement noise ---
    sigma_meas = 5 * ones(4,1);  % [cm]

    % --- Actuator limits ---
    u_min = [0; 0];             % [cm³/s]
    u_max = [1000; 1000];       % [cm³/s]

    % --- Disturbance setup (piecewise constant every 60s) ---
    Delta_t_d = 60;                  
    t_d = 0:Delta_t_d:t(end);
    N_d = numel(t_d);
    rng(2);  % reproducibility
    d_values = 200 + 100*rand(2, N_d);   % each block ∈ [200, 300] cm³/s

    % Zero-order hold interpolants
    g1 = griddedInterpolant(t_d, d_values(1,:), 'previous');
    g2 = griddedInterpolant(t_d, d_values(2,:), 'previous');
    d_of_t = @(tt) [g1(tt); g2(tt)];

    % ------------------------------------------------------------
    % Simulation loop (Euler integration + PID)
    % ------------------------------------------------------------
    for k = 1:N-1
        % Current disturbance
        d_current = d_of_t(t(k));

        % Sensor readings
        y_true(:,k) = FourTankSystemSensor(x(:,k), p);
        y_meas(:,k) = y_true(:,k) + randn(4,1) .* sigma_meas;

        % Current reference (for tanks 1 and 2)
        h_ref_current = h_ref_traj(:,k);

        % --- PID control law (deviation variables) ---
        e = h_ref_current - y_meas(1:2,k);
        int_e = int_e + e * Ts;
        der_e = (e - prev_e) / Ts;
        prev_e = e;

        du = Kp .* e + Ki .* int_e + Kd .* der_e;   % control deviation
        u(:,k) = u_s + du;                          % absolute control

        % --- Enforce actuator limits ---
        u(:,k) = max(min(u(:,k), u_max), u_min);

        % --- System update (Euler) ---
        f = FourTankSystemModified(0, x(:,k), u(:,k), d_current, p);
        x(:,k+1) = x(:,k) + f * Ts;
    end

    % --- Final outputs ---
    y_true(:,N) = FourTankSystemSensor(x(:,N), p);
    y_meas(:,N) = y_true(:,N) + randn(4,1).*sigma_meas;

    % ------------------------------------------------------------
    % Plot results
    % ------------------------------------------------------------
        % --- PID string for figure titles ---
    pid_str = sprintf('Kp=[%.2f,%.2f], Ki=[%.2f,%.2f], Kd=[%.2f,%.2f]', ...
                       Kp(1),Kp(2),Ki(1),Ki(2),Kd(1),Kd(2));

    % Heights vs references
    [fig1,~] = plot_heights_2(y_meas, y_true, t, p, h_ref_traj);
    sgtitle(fig1, 'Tank Heights – PID Control (Stochastic, piecewise disturbances)');
    exportgraphics(fig1, fullfile(outputFolder, 'heights_pid_stochastic.pdf'), 'ContentType', 'vector');

    % Inputs with max line
    fig2 = figure;
    hold on; grid on;
    plot(t, u(1,:), 'LineWidth', 1.5, 'DisplayName', 'Pump 1');
    plot(t, u(2,:), 'LineWidth', 1.5, 'DisplayName', 'Pump 2');
    yline(u_max(1), '--r', 'Max Pump Flow (1000 cm³/s)', 'LineWidth', 1.2);
    xlabel('Time [s]');
    ylabel('Flow [cm³/s]');
    title('Control Inputs (PID – Stochastic, with saturation)');
    legend('Location','best');
    ylim([0 u_max(1)*1.1]);
    xlim([t(1) t(end)])
    exportgraphics(fig2, fullfile(outputFolder, 'inputs_pid_stochastic.pdf'), 'ContentType', 'vector');

    % Disturbance realization
    fig3 = figure;
    hold on; grid on;
    stairs(t_d, d_values(1,:), 'LineWidth', 1.5, 'DisplayName', 'F_3');
    stairs(t_d, d_values(2,:), 'LineWidth', 1.5, 'DisplayName', 'F_4');
    xlabel('Time [s]');
    ylabel('Disturbance Flow [cm³/s]');
    title('Piecewise Constant Disturbances');
    legend('Location','best');
    xlim([t(1), t(end)]);
    exportgraphics(fig3, fullfile(outputFolder, 'disturbances_pid_stochastic.pdf'), 'ContentType', 'vector');

    % 4️⃣ Inputs and Disturbances (stacked figure)
    fig = figure('Name','Inputs and Disturbances');
    % --- Upper subplot: Inputs ---
    subplot(2,1,1);
    hold on; grid on;
    stairs(t, u(1,:), 'LineWidth', 1.5, 'DisplayName', 'F_1 (Pump 1)');
    stairs(t, u(2,:), 'LineWidth', 1.5, 'DisplayName', 'F_2 (Pump 2)');
    yline(u_max(1), '--r', 'DisplayName', 'Flow limit (1000 cm³/s)', 'LineWidth', 1.2);
    ylabel('Input Flow [cm³/s]');
    title(sprintf('Control Inputs (PID – Stochastic)\n%s', pid_str));
    legend('Location','best');
    xlim([t(1), t(end)]);
    ylim([0, 1000])

    % --- Lower subplot: Disturbances ---
    subplot(2,1,2);
    hold on; grid on;
    stairs(t_d, d_values(1,:), 'LineWidth', 1.5, 'DisplayName', 'F_3');
    stairs(t_d, d_values(2,:), 'LineWidth', 1.5, 'DisplayName', 'F_4');
    xlabel('Time [s]');
    ylabel('Disturbance Flow [cm³/s]');
    title('Piecewise Constant Disturbances (F₃ and F₄)');
    legend('Interpreter','latex','Location','best');
    xlim([t(1), t(end)]);
    ylim([0, 500]);

    %sgtitle('Inputs and Disturbances (Stochastic Simulation)');
    %sgtitle(sprintf('Control Inputs (PID – Stochastic)\n%s', pid_str));
    exportgraphics(fig, fullfile(outputFolder, 'inputs_and_disturbances.pdf'), 'ContentType', 'vector');

    disp(['✅ PID stochastic simulation complete. Results saved in: ', outputFolder]);
end
