function [fig, ax] = plot_step_responses(t, y_10, y_25, y_50, p, stepSource)
%PLOT_STEP_RESPONSES Plot nonlinear step responses for all 4 tanks
%
%   [fig, ax] = plot_step_responses(t, y_10, y_25, y_50, p, stepSource)
%
%   Inputs:
%       t           : time vector [s]
%       y_10, y_25, y_50 : tank heights or masses for 10%, 25%, 50% steps (4xN)
%       p           : parameter vector (unused, for consistency)
%       stepSource  : string, e.g. 'F_1' or 'F_2'
%
%   Output:
%       fig, ax : figure and axes handles
%
%   Example:
%       plot_step_responses(t, y10, y25, y50, p, 'F_1');

    if nargin < 6
        stepSource = 'F_1'; % default if not specified
    end

    % --- Plot layout ---
    order = [3, 4, 1, 2];   % consistent tank order
    h_max = max([y_10(:); y_25(:); y_50(:)]);
    h_min = min([y_10(:); y_25(:); y_50(:)]);

    % --- Create figure ---
    fig = figure('Name', sprintf('Step Responses for %s', stepSource));
    sgtitle(sprintf('Tank Responses to Step Changes in %s', stepSource), ...
            'Interpreter','latex', 'FontSize',14, 'FontWeight','bold');

    ax = gobjects(1,4);
    colors = lines(3); % for consistent coloring

    for i = 1:4
        ax(i) = subplot(2,2,order(i)); hold on; grid on;

        % Plot step responses
        plot(t, y_10(i,:), 'Color', colors(1,:), 'LineWidth', 1.5, ...
             'DisplayName', '10\% Step');
        plot(t, y_25(i,:), 'Color', colors(2,:), 'LineWidth', 1.5, ...
             'DisplayName', '25\% Step');
        plot(t, y_50(i,:), 'Color', colors(3,:), 'LineWidth', 1.5, ...
             'DisplayName', '50\% Step');

        xlabel('Time [s]');
        ylabel(sprintf('$h_{%d}$ [cm]', i), 'Interpreter', 'latex');
        title(sprintf('Tank %d', i), 'FontWeight','bold');
        legend('Location', 'best', 'Interpreter','latex');
        xlim([t(1), t(end)]);
        ylim([0.9*h_min 1.1*h_max]);
    end

    % Adjust layout
    set(gcf, 'Position', [100 100 900 700]);
end
