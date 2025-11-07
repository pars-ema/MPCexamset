function [fig, ax] = plot_transfer_identification_combined(t, y_norm_F1_10, y_norm_F1_25, y_norm_F1_50, ...
                                                           y_norm_F2_10, y_norm_F2_25, y_norm_F2_50, p)
%PLOT_TRANSFER_IDENTIFICATION_COMBINED Plot normalized transfer identification for multiple step sizes
%
%   [fig, ax] = plot_transfer_identification_combined(t, y_norm_F1_10, y_norm_F1_25, y_norm_F1_50, ...
%                                                     y_norm_F2_10, y_norm_F2_25, y_norm_F2_50, p)
%
%   Inputs:
%       t                : time vector [s]
%       y_norm_F1_*      : normalized responses for F1 step (4×N)
%       y_norm_F2_*      : normalized responses for F2 step (4×N)
%       p                : parameter vector (optional)
%
%   Output:
%       fig, ax          : figure and axes handles (2×2 grid)
%
%   Each subplot shows a transfer G_ij(s):
%       G11: h1 ← F1 | G12: h1 ← F2
%       G21: h2 ← F1 | G22: h2 ← F2

    % --- Setup ---
    fig = figure('Name', 'Normalized Step Responses for Transfer Identification (All Steps)');
    sgtitle({'Normalized Step Responses'}, ...
             'Interpreter','latex', 'FontSize', 15, 'FontWeight','bold');

    colors = lines(3);
    ax = gobjects(2,2);
    orderNames = {'$G_{11}: h_1 \leftarrow F_1$', '$G_{12}: h_1 \leftarrow F_2$', ...
                  '$G_{21}: h_2 \leftarrow F_1$', '$G_{22}: h_2 \leftarrow F_2$'};

    % --- Determine global axis limits ---
    allData = [y_norm_F1_10(:); y_norm_F1_25(:); y_norm_F1_50(:); ...
               y_norm_F2_10(:); y_norm_F2_25(:); y_norm_F2_50(:)];
    y_min = min(allData);
    y_max = max(allData);

    % --- Subplots ---
    ax(1,1) = subplot(2,2,1); hold on; grid on;
    plot(t, y_norm_F1_10(1,:), 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', '10\% Step');
    plot(t, y_norm_F1_25(1,:), 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', '25\% Step');
    plot(t, y_norm_F1_50(1,:), 'Color', colors(3,:), 'LineWidth', 2, 'DisplayName', '50\% Step');
    title(orderNames{1}, 'Interpreter','latex', 'FontSize',14,'FontWeight','bold');

    ax(1,2) = subplot(2,2,2); hold on; grid on;
    plot(t, y_norm_F2_10(1,:), 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', '10\% Step');
    plot(t, y_norm_F2_25(1,:), 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', '25\% Step');
    plot(t, y_norm_F2_50(1,:), 'Color', colors(3,:), 'LineWidth', 2, 'DisplayName', '50\% Step');
    title(orderNames{2}, 'Interpreter','latex', 'FontSize',14,'FontWeight','bold');

    ax(2,1) = subplot(2,2,3); hold on; grid on;
    plot(t, y_norm_F1_10(2,:), 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', '10\% Step');
    plot(t, y_norm_F1_25(2,:), 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', '25\% Step');
    plot(t, y_norm_F1_50(2,:), 'Color', colors(3,:), 'LineWidth', 2, 'DisplayName', '50\% Step');
    title(orderNames{3}, 'Interpreter','latex', 'FontSize',14,'FontWeight','bold');

    ax(2,2) = subplot(2,2,4); hold on; grid on;
    plot(t, y_norm_F2_10(2,:), 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', '10\% Step');
    plot(t, y_norm_F2_25(2,:), 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', '25\% Step');
    plot(t, y_norm_F2_50(2,:), 'Color', colors(3,:), 'LineWidth', 2, 'DisplayName', '50\% Step');
    title(orderNames{4}, 'Interpreter','latex', 'FontSize',14,'FontWeight','bold');

    % --- Axes styling ---
    for i = 1:numel(ax)
        xlabel(ax(i), 'Time [s]', 'FontSize', 14, 'FontWeight','bold');
        ylabel(ax(i), 'Normalized Response', 'Interpreter','latex', 'FontSize', 14, 'FontWeight','bold');
        %ylim(ax(i), [y_min*1.1, y_max*1.1]);
        ylim(ax(i), [-0.2, 0.6]);
        xlim(ax(i), [t(1), t(end)]);
        ax(i).FontSize = 12;
        ax(i).LineWidth = 1.2;
        legend(ax(i), 'Location','best', 'Interpreter','latex', 'FontSize', 12, 'Box','off');
    end

    % --- Adjust layout ---
    set(fig, 'Position', [100 100 950 750]);
end
