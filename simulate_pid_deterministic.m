function simulate_pid_deterministic(t, x0, u_s, d, p, h_ref_traj, Kp, Ki, Kd, outputFolder)
% SIMULATE_PID_DETERMINISTIC
% Closed-loop PID control for the nonlinear four-tank system (deterministic model)
% with actuator limits and PID values shown in figure titles.

    if nargin < 10 || isempty(outputFolder)
        outputFolder = fullfile('figures','problem_3','pid_deterministic');
    end
    if ~exist(outputFolder,'dir'); mkdir(outputFolder); end

    Ts = t(2) - t(1);
    N  = length(t);
    nx = 4; ny = 4; nu = 2;
    x = zeros(nx,N);
    y_true = zeros(ny,N);
    y_meas = zeros(ny,N);
    u = zeros(nu,N);
    int_e = zeros(2,1);
    prev_e = zeros(2,1);
    x(:,1) = x0;

    sigma_meas = zeros(4,1); % deterministic → no noise
    u_min = [0;0]; 
    u_max = [1000;1000]; % [cm³/s]

    % --- Simulation loop ---
    for k = 1:N-1
        y_true(:,k) = FourTankSystemSensor(x(:,k), p);
        y_meas(:,k) = y_true(:,k) + randn(4,1).*sigma_meas;

        h_ref_current = h_ref_traj(:,k);
        e = h_ref_current - y_meas(1:2,k);
        int_e = int_e + e * Ts;
        der_e = (e - prev_e) / Ts;
        prev_e = e;

        du = Kp .* e + Ki .* int_e + Kd .* der_e;
        u(:,k) = u_s + du;
        u(:,k) = max(min(u(:,k), u_max), u_min);

        f = FourTankSystemModified(0, x(:,k), u(:,k), d, p);
        x(:,k+1) = x(:,k) + f * Ts;
    end

    y_true(:,N) = FourTankSystemSensor(x(:,N), p);
    y_meas(:,N) = y_true(:,N);

    % --- PID string for figure titles ---
    pid_str = sprintf('Kp=[%.2f,%.2f], Ki=[%.2f,%.2f], Kd=[%.2f,%.2f]', ...
                       Kp(1),Kp(2),Ki(1),Ki(2),Kd(1),Kd(2));

    % ------------------------------------------------------------
    % Plot results
    % ------------------------------------------------------------
    [fig1,~] = plot_heights_2(y_meas, y_true, t, p, h_ref_traj);
    sgtitle(fig1, sprintf('Tank Heights – PID Control (Deterministic)\n%s', pid_str));
    exportgraphics(fig1, fullfile(outputFolder, 'heights_pid_deterministic.pdf'), 'ContentType', 'vector');

    fig2 = figure;
    hold on; grid on;
    plot(t, u(1,:), 'LineWidth', 1.5, 'DisplayName', 'Pump 1');
    plot(t, u(2,:), 'LineWidth', 1.5, 'DisplayName', 'Pump 2');
    yline(u_max(1), '--r', 'Max Pump Flow (1000 cm³/s)', 'LineWidth', 1.2);
    xlabel('Time [s]');
    ylabel('Flow [cm³/s]');
    title(sprintf('Control Inputs (PID – Deterministic)\n%s', pid_str));
    legend('Location','best');
    ylim([0 u_max(1)*1.1]);
    xlim([t(1) t(end)])
    exportgraphics(fig2, fullfile(outputFolder, 'inputs_pid_deterministic.pdf'), 'ContentType', 'vector');

        % 3️⃣ Inputs and Disturbances (stacked vertically)
    fig = figure('Name','Inputs and Disturbances');
    
    % --- Upper subplot: Inputs ---
    subplot(2,1,1);
    hold on; grid on;
    stairs(t, u(1,:), 'LineWidth', 1.5, 'DisplayName', 'F_1 (Pump 1)');
    stairs(t, u(2,:), 'LineWidth', 1.5, 'DisplayName', 'F_2 (Pump 2)');
    yline(u_max(1), '--r', 'DisplayName', 'Flow limit (1000 cm³/s)', ...
      'LineWidth', 1.2);
    ylabel('Input Flow [cm³/s]');
    title(sprintf('Control Inputs (PID – Deterministic)\n%s', pid_str));
    legend('Interpreter','latex','Location','best');
    xlim([t(1), t(end)]);
    ylim([0 u_max(1)*1.1]);
    
    % --- Lower subplot: Disturbances ---
    d = repmat(d, length(t));
    d;
    subplot(2,1,2);
    hold on; grid on;
    stairs(t, d(1,:), 'LineWidth', 1.5, 'DisplayName', 'F_3');
    stairs(t, d(2,:), 'LineWidth', 1.5, 'DisplayName', 'F_4');
    xlabel('Time [s]');
    ylabel('Disturbance Flow [cm³/s]');
    title('Disturbances (Unmeasured F₃ and F₄)');
    legend('Location','best');
    xlim([t(1), t(end)]);

    %sgtitle('Inputs and Disturbances (Deterministic Simulation)');
    %sgtitle(sprintf('Control Inputs (PID – Deterministic)\n%s', pid_str), fontsize=11);



    exportgraphics(fig, fullfile(outputFolder, 'inputs_and_disturbances.pdf'), 'ContentType', 'vector');

    disp(['✅ PID deterministic simulation complete. Results saved in: ', outputFolder]);
end
