function [x, y_true_abs, y_meas_abs, y_meas_dev] = simulate_sde_2( ...
        t, x0, u, d, p, sigma_meas, show_plot, outputFolder, output_index, xs)

% SIMULATE_SDE_2
% Stochastic nonlinear Four-Tank System (Euler–Maruyama) with
% deviation measurements for Kalman filtering.
%
% Inputs:
%   t           : time vector [1×N]
%   x0          : initial state [4×1]
%   u           : inputs [2×1] or [2×N]
%   d           : disturbances [2×1] or [2×N]
%   p           : parameter struct/vector
%   sigma_meas  : measurement noise std [4×1]
%   show_plot   : logical, plot results or not
%   outputFolder: folder for plots (optional)
%   output_index: measured output indices (e.g. [1 2])
%   xs          : steady-state state [4×1] (for deviation)
%
% Outputs:
%   x            : states [4×N]
%   y_true_abs   : true outputs (all 4 tanks) [4×N]
%   y_meas_abs   : measured outputs (subset) [ny_meas×N]
%   y_meas_dev   : deviation measurements (subset) [ny_meas×N]

    % ----------------------------------------------------
    % Defaults
    % ----------------------------------------------------
    if nargin < 7 || isempty(show_plot), show_plot = false; end
    if nargin < 8 || isempty(outputFolder), outputFolder = fullfile('figures','sde'); end
    if nargin < 9 || isempty(output_index), output_index = 1:4; end
    if nargin < 10 || isempty(xs), xs = x0; end

    Ts = t(2) - t(1);
    N  = length(t);
    nx = 4; ny = 4;
    ny_meas = length(output_index);

    % Preallocate
    x           = zeros(nx,N);
    y_true_abs  = zeros(ny,N);
    y_meas_abs  = zeros(ny_meas,N);
    y_meas_dev  = zeros(ny_meas,N);

    x(:,1) = x0;

    % Make u and d 2×N
    if size(u,2) == 1, u = repmat(u,1,N); end
    if size(d,2) == 1, d = repmat(d,1,N); end

    % Wiener noise for F3,F4
    tf = t(end);
    [~,~,dW3] = ScalarStdWienerProcess(tf, N, 1);
    [~,~,dW4] = ScalarStdWienerProcess(tf, N, 1);
    sigma3 = 20;
    sigma4 = 20;

    % Steady-state output (absolute) and its measured part
    y_s_abs = FourTankSystemSensor(xs,p);
    y_s_meas = y_s_abs(output_index);

    % Initial output
    y_true_abs(:,1) = FourTankSystemSensor(x(:,1),p);
    y_meas_abs(:,1) = y_true_abs(output_index,1) + ...
                      sigma_meas(output_index).*randn(ny_meas,1);
    y_meas_dev(:,1) = y_meas_abs(:,1) - y_s_meas;

    % ----------------------------------------------------
    % Euler–Maruyama loop
    % ----------------------------------------------------
    for k = 1:N-1

        f = FourTankSystemModified(0, x(:,k), u(:,k), d(:,k), p);
        noise = [0;
                 0;
                 sigma3*dW3(k);
                 sigma4*dW4(k)];

        x(:,k+1) = x(:,k) + Ts*f + noise;

        y_true_abs(:,k+1) = FourTankSystemSensor(x(:,k+1),p);
        y_meas_abs(:,k+1) = y_true_abs(output_index,k+1) + ...
                            sigma_meas(output_index).*randn(ny_meas,1);
        y_meas_dev(:,k+1) = y_meas_abs(:,k+1) - y_s_meas;
    end

    % ----------------------------------------------------
    % Plotting (absolute + measured)
    % ----------------------------------------------------
    if show_plot
        if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
        T = t;

        % Masses
        [fig1,~] = plot_masses(x,T,p);
        sgtitle(fig1,'Masses (SDE Model)');
        exportgraphics(fig1, fullfile(outputFolder,'masses_SDE.pdf'), 'ContentType','vector');

        % Heights
        [fig2,~] = plot_heights(x,T,p);
        sgtitle(fig2,'Heights (SDE Model)');
        exportgraphics(fig2, fullfile(outputFolder,'heights_SDE.pdf'), 'ContentType','vector');

        % Measured vs True (absolute)
        Ypad = pad_measured(y_meas_abs, output_index, N);
        [fig3,~] = plot_heights_2(Ypad, y_true_abs, T, p);
        sgtitle(fig3,'Measured vs True Heights (SDE Model)');
        exportgraphics(fig3, fullfile(outputFolder,'heights_measured_vs_true_SDE.pdf'), 'ContentType','vector');
    end
end

%% helper: expand measured subset to 4×N for plotting only
function Ypad = pad_measured(y_meas_sel, idx, N)
    Ypad = nan(4,N);
    Ypad(idx,:) = y_meas_sel;
end
