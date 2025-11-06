function simulate_stochastic(t, x0, u, p, outputFolder)
% SIMULATE_STOCHASTIC
% Simulate stochastic nonlinear model with piecewise-constant disturbances.
% Disturbances F3, F4 vary every 60 s, measurement noise added.
%
% Inputs:
%   t  - time vector [s]
%   x0 - initial state [4x1]
%   u  - input sequence [2xN] or constant [2x1]
%   p  - parameter struct
%   outputFolder - folder path for saving figures (optional)
% -------------------------------------------------------------

    % --- Output folder ---
    if nargin < 5 || isempty(outputFolder)
        outputFolder = fullfile('figures', 'stochastic');
    end
    if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end

    % --- Disturbance setup ---
    Delta_t_d = 60;              % [s] disturbance change interval
    t_d = 0:Delta_t_d:t(end);
    N_d = numel(t_d);

    rng(1);                      % reproducible seed
    d_values = 200 + 100*rand(2, N_d);  % each block in [200, 300] cm³/s

    g1 = griddedInterpolant(t_d, d_values(1,:), 'previous');
    g2 = griddedInterpolant(t_d, d_values(2,:), 'previous');
    d_of_t = @(tt) [g1(tt); g2(tt)];    % returns 2x1 disturbance vector

    % --- Time setup ---
    Ts = t(2) - t(1);
    N  = numel(t);

    % --- State allocation ---
    nx = 4; ny = 4;
    x = zeros(nx, N);
    y_true = zeros(ny, N);
    y_meas = zeros(ny, N);
    sigma_y = 5 * ones(ny, 1);  % measurement noise [cm]
    x(:,1) = x0;

    % --- Handle constant vs time-varying input ---
    if size(u,2) == 1
        u = repmat(u, 1, N);
    end

    % --- Euler integration loop ---
    for k = 1:N-1
        d_current = d_of_t(t(k));
        f = FourTankSystemModified(0, x(:,k), u(:,k), d_current, p);
        x(:,k+1) = x(:,k) + f * Ts;

        % Sensor model
        y_true(:,k+1) = FourTankSystemSensor(x(:,k+1), p);
        y_meas(:,k+1) = y_true(:,k+1) + sigma_y .* randn(4,1);
    end

    % --- Prepare data for plotting ---
    X = x;
    T = t;

    % ------------------------------------------------------------
    % Plot results
    % ------------------------------------------------------------
    % 1️⃣ Heights (true vs measured)
    [fig3, ~] = plot_heights_2(y_meas, y_true, T, p);
    sgtitle(fig3, 'Measured vs True Heights (Stochastic Nonlinear Model)');
    exportgraphics(fig3, fullfile(outputFolder, 'heights_measured_vs_true.pdf'), 'ContentType', 'vector');

    % 2️⃣ Tank heights (true)
    [fig2, ~] = plot_heights(X, T, p);
    sgtitle(fig2, 'Heights in the tanks (Stochastic Nonlinear Model)');
    exportgraphics(fig2, fullfile(outputFolder, 'heights_stochastic.pdf'), 'ContentType', 'vector');

    % 3️⃣ Masses
    [fig1, ~] = plot_masses(X, T, p);
    sgtitle(fig1, 'Masses in the tanks (Stochastic Nonlinear Model)');
    exportgraphics(fig1, fullfile(outputFolder, 'masses_stochastic.pdf'), 'ContentType', 'vector');

    % 4️⃣ Inputs and Disturbances (stacked figure)
    fig = figure('Name','Inputs and Disturbances');
    % --- Upper subplot: Inputs ---
    subplot(2,1,1);
    hold on; grid on;
    stairs(T, u(1,:), 'LineWidth', 1.5, 'DisplayName', 'F_1 (Pump 1)');
    stairs(T, u(2,:), 'LineWidth', 1.5, 'DisplayName', 'F_2 (Pump 2)');
    ylabel('Input Flow [cm³/s]');
    title('Inputs (Pumps F₁ and F₂)');
    legend('Location','best');
    xlim([T(1), T(end)]);
    ylim([0, 1000])

    % --- Lower subplot: Disturbances ---
    subplot(2,1,2);
    hold on; grid on;
    stairs(t_d, d_values(1,:), 'LineWidth', 1.5, 'DisplayName', 'F_3');
    stairs(t_d, d_values(2,:), 'LineWidth', 1.5, 'DisplayName', 'F_4');
    xlabel('Time [s]');
    ylabel('Disturbance Flow [cm³/s]');
    title('Piecewise Constant Disturbances (F₃ and F₄)');
    legend('Location','best');
    xlim([T(1), T(end)]);
    ylim([0, 500]);

    sgtitle('Inputs and Disturbances (Stochastic Simulation)');
    exportgraphics(fig, fullfile(outputFolder, 'inputs_and_disturbances.pdf'), 'ContentType', 'vector');

    disp(['✅ Stochastic simulation complete. Results saved in: ', outputFolder]);
end
