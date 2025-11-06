function simulate_sde(t, x0, u, d, p, outputFolder)
% SIMULATE_SDE
% Stochastic nonlinear model using Euler–Maruyama integration.
% Matches the behavior of the working SDE script version.
%
% Inputs:
%   t  - time vector [s]
%   x0 - initial state [4x1]
%   u  - input sequence [2xN] or constant [2x1]
%   d  - mean disturbance [2x1]
%   p  - parameter vector or struct
%   outputFolder - path for saving figures
%
% Author: MPC Project Team
% ------------------------------------------------------------

    if nargin < 6 || isempty(outputFolder)
        outputFolder = fullfile('figures','sde');
    end
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end

    % --- Time setup ---
    Ts = t(2) - t(1);
    N  = length(t);
    Delta_t = Ts;

    % --- Initialization ---
    nx = 4; ny = 4; nz = 2;
    x = zeros(nx, N);
    y_true = zeros(ny, N);
    y_measured = zeros(ny, N);
    z = zeros(nz, N);
    T = zeros(1, N);
    x(:,1) = x0;

    % --- Generate Wiener processes ---
    tf = t(end);
    [W3, TW3, dW3] = ScalarStdWienerProcess(tf, N, 1); % for F3
    [W4, TW4, dW4] = ScalarStdWienerProcess(tf, N, 1); % for F4

    % --- Noise intensities ---
    sigma1 = 0;
    sigma2 = 0;
    sigma3 = 100;
    sigma4 = 100;

    % --- Mean disturbances ---
    F3 = d(1);
    F4 = d(2);

    % --- Plot Wiener processes and noisy disturbances ---
    figW = figure;
    hold on; grid on;
    plot(TW3, W3, 'DisplayName', 'W3');
    plot(TW4, W4, 'DisplayName', 'W4');
    plot(TW3(2:end), F3 + sigma3*dW3, 'DisplayName', 'F3 + σ₃·dW₃');
    plot(TW4(2:end), F4 + sigma4*dW4, 'DisplayName', 'F4 + σ₄·dW₄');
    yline(F3, '--', 'DisplayName', 'Mean F3');
    yline(F4, '--', 'DisplayName', 'Mean F4');
    xlabel('Time [s]');
    ylabel('Wiener process value');
    title('Two Independent Wiener Processes');
    legend show;
    exportgraphics(figW, fullfile(outputFolder, 'WD_SDE.pdf'), 'ContentType', 'vector');
    ylim([0, 1000])

    % --- Measurement noise standard deviations ---
    sigma_measurement = [5; 5; 5; 5];

    % --- Simulation (Euler–Maruyama) ---
    for k = 1:N-1
        % Deterministic dynamics
        f = FourTankSystemModified(0, x(:,k), u(:,min(k, size(u,2))), d, p);

        % Stochastic term (affecting tank 3 and 4)
        noise = [0; 0; sigma3*dW3(k); sigma4*dW4(k)];
        % Euler–Maruyama update
        x(:,k+1) = x(:,k) + f * Delta_t + noise;
        % Outputs
        y_true(:,k+1) = FourTankSystemSensor(x(:,k+1), p);
        y_measured(:,k+1) = y_true(:,k+1) + randn(4,1) .* sigma_measurement;
        z(:,k+1) = FourTankSystemOutput(x(:,k+1), p);
        T(k+1) = (k+1) * Delta_t;
    end

    % ------------------------------------------------------------
    % Plot results
    % ------------------------------------------------------------
    % 1️⃣ Masses
    [fig1, ~] = plot_masses(x, T, p);
    sgtitle(fig1, 'Masses in the tanks (SDE Model)');
    exportgraphics(fig1, fullfile(outputFolder, 'masses_SDE.pdf'), 'ContentType', 'vector');

    % 2️⃣ Heights
    [fig2, ~] = plot_heights(x, T, p);
    sgtitle(fig2, 'Heights in the tanks (SDE Model)');
    exportgraphics(fig2, fullfile(outputFolder, 'heights_SDE.pdf'), 'ContentType', 'vector');

    % 3️⃣ Measured vs True Heights
    [fig3, ~] = plot_heights_2(y_measured, y_true, T, p);
    sgtitle(fig3, 'Measured vs True Heights (SDE Model)');
    exportgraphics(fig3, fullfile(outputFolder, 'heights_measured_vs_true_SDE.pdf'), 'ContentType', 'vector');

    % 4️⃣ Inputs + Disturbances (stacked)
    fig4 = figure('Name','Inputs and Disturbances');
    subplot(2,1,1);
    hold on; grid on;
    if size(u,2) == 1
        stairs(T, u(1)*ones(size(T)), 'LineWidth', 1.5, 'DisplayName', 'F_1');
        stairs(T, u(2)*ones(size(T)), 'LineWidth', 1.5, 'DisplayName', 'F_2');
    else
        stairs(T, u(1,:), 'LineWidth', 1.5, 'DisplayName', 'F_1');
        stairs(T, u(2,:), 'LineWidth', 1.5, 'DisplayName', 'F_2');
    end
    ylabel('Input Flow [cm³/s]');
    title('Inputs (Pumps F₁ and F₂)');
    legend('Location','best');
    xlim([T(1), T(end)]);

    subplot(2,1,2);
    hold on; grid on;
    plot(TW3(2:end), F3 + sigma3*dW3, 'LineWidth', 1.5, 'DisplayName', 'F₃ + σ₃·dW₃');
    plot(TW4(2:end), F4 + sigma4*dW4, 'LineWidth', 1.5, 'DisplayName', 'F₄ + σ₄·dW₄');
    yline(F3, '--', 'Mean F₃');
    yline(F4, '--', 'Mean F₄');
    xlabel('Time [s]');
    ylabel('Disturbance Flow [cm³/s]');
    title('Disturbances (F₃, F₄ with SDE noise)');
    legend('Location','best');
    xlim([T(1), T(end)]);
    ylim([0, 1000])

    sgtitle('Inputs and Disturbances (SDE Simulation)');
    exportgraphics(fig4, fullfile(outputFolder, 'inputs_and_disturbances_SDE.pdf'), 'ContentType', 'vector');

    disp(['✅ SDE simulation complete. Results saved in: ', outputFolder]);
end
