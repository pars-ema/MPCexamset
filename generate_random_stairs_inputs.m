function [u1, u2, t_change_1, t_change_2, values_1, values_2] = ...
    generate_random_stairs_inputs(t, n_changes, Fmin, Fmax, seed)
% GENERATE_RANDOM_STAIRS_INPUTS
% ---------------------------------------------------------------
% Generates two random piecewise-constant (staircase) input signals.
%
% Inputs:
%   t          - time vector [s]
%   n_changes  - number of step changes for each input
%   Fmin       - minimum flow value [cm³/s]
%   Fmax       - maximum flow value [cm³/s]
%   seed       - RNG seed for reproducibility
%
% Outputs:
%   u1, u2         - input signals (vectors same size as t)
%   t_change_1,2   - time instants of value changes for each input
%   values_1,2     - corresponding flow values for each segment
%
% Example:
%   t = 0:1:3600;
%   [u1,u2,tc1,tc2,v1,v2] = generate_random_stairs_inputs(t,5,100,300,1);
%
% Author: ChatGPT (refactor for MPC project)
% ---------------------------------------------------------------

    if nargin < 5
        seed = 1; % default seed
    end
    rng(seed);

    t0 = t(1);
    tf = t(end);
    Nt = length(t);

    % --- Random change instants (sorted unique indices) ---
    t_change_1 = sort(randperm(Nt-2, n_changes));
    t_change_1 = t(t_change_1);
    t_change_1 = [t0, t_change_1, tf]; % include start and end

    t_change_2 = sort(randperm(Nt-2, n_changes));
    t_change_2 = t(t_change_2);
    t_change_2 = [t0, t_change_2, tf];

    % --- Random values (within given range) ---
    values_1 = Fmin + (Fmax - Fmin) * rand(1, length(t_change_1));
    values_2 = Fmin + (Fmax - Fmin) * rand(1, length(t_change_2));

    % --- Use your existing helper to make step signals ---
    [u1, u2] = generate_stairs_inputs(t, t_change_1, values_1, t_change_2, values_2);

    % --- Optional plot for quick check ---
    figure;
    stairs(t, u1, 'LineWidth', 1.5, 'DisplayName','Pump 1');
    hold on;
    stairs(t, u2, 'LineWidth', 1.5, 'DisplayName','Pump 2');
    xlabel('Time [s]');
    ylabel('Flow [cm³/s]');
    title(sprintf('Random Staircase Inputs (%d changes each)', n_changes));
    legend('Location','best');
    grid on;
    xlim([t0, tf]);

end
