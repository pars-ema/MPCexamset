function [t, y] = simulate_deterministic(t, x0, u, d, p, plot, outputFolder)
% SIMULATE_DETERMINISTIC
% Simulate deterministic nonlinear model using a fixed-step Euler integrator.
% Also plots u (inputs) and d (disturbances) in one figure, stacked vertically.
%
% Inputs:
%   t  - time vector [s]
%   x0 - initial state [4x1]
%   u  - input sequence [2xN] or constant [2x1]
%   d  - disturbance sequence [2xN] or constant [2x1]
%   p  - parameter struct
%   outputFolder - folder path for saving figures
% -------------------------------------------------------------

    % --- Output folder ---
    if plot
    if nargin < 6 || isempty(outputFolder)
        outputFolder = fullfile('figures', 'deterministic');
    end
    if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
    end

    % --- Time setup ---
    Ts = t(2) - t(1);
    N  = numel(t);

    % --- State allocation ---
    nx = 4;
    x = zeros(nx, N);
    x(:,1) = x0;

    % --- Handle constant vs time-varying input/disturbance ---
    if size(u,2) == 1
        u = repmat(u, 1, N);
    end
    if size(d,2) == 1
        d = repmat(d, 1, N);
    end

    % --- Euler integration loop ---
    for k = 1:N-1
        f = FourTankSystemModified(t(k), x(:,k), u(:,k), d(:,k), p);
        x(:,k+1) = x(:,k) + f * Ts;
    end

    % --- Convert for plotting ---
    X = x;
    T = t;
    y = FourTankSystemSensor(x, p)
    if plot
        % ------------------------------------------------------------
        % Plot results
        % ------------------------------------------------------------
        % 1️⃣ Masses
        [fig, ~] = plot_masses(X, T, p);
        sgtitle(fig, 'Masses in the tanks (Deterministic Model)');
        exportgraphics(fig, fullfile(outputFolder, 'masses.pdf'), 'ContentType', 'vector');
        
        % 2️⃣ Heights
        [fig, ~] = plot_heights(X, T, p);
        sgtitle(fig, 'Heights in the tanks (Deterministic Model)');
        exportgraphics(fig, fullfile(outputFolder, 'heights.pdf'), 'ContentType', 'vector');
        
        % 3️⃣ Inputs and Disturbances (stacked vertically)
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
        stairs(T, d(1,:), 'LineWidth', 1.5, 'DisplayName', 'F_3');
        stairs(T, d(2,:), 'LineWidth', 1.5, 'DisplayName', 'F_4');
        xlabel('Time [s]');
        ylabel('Disturbance Flow [cm³/s]');
        title('Disturbances (Unmeasured F₃ and F₄)');
        legend('Location','best');
        xlim([T(1), T(end)]);
        
        sgtitle('Inputs and Disturbances (Deterministic Simulation)');
        set(gcf,'PaperPositionMode','auto');
        
        exportgraphics(fig, fullfile(outputFolder, 'inputs_and_disturbances.pdf'), 'ContentType', 'vector');
        
        disp(['✅ Deterministic simulation complete. Results saved in: ', outputFolder]);
    end
end
