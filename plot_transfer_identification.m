function [fig, ax] = plot_transfer_identification(t, y_norm_F1, y_norm_F2)
%PLOT_TRANSFER_IDENTIFICATION  Plot normalized step responses in 2x2 grid
%
%   [fig, ax] = plot_transfer_identification(t, y_norm_F1, y_norm_F2)
%
%   Inputs:
%       t          : time vector [s]
%       y_norm_F1  : normalized step responses for F1 step (4×N)
%       y_norm_F2  : normalized step responses for F2 step (4×N)
%
%   Outputs:
%       fig, ax : figure handle and axes handles (2×2)
%
%   Layout:
%       G11 : h1 ← F1   |  G12 : h1 ← F2
%       G21 : h2 ← F1   |  G22 : h2 ← F2

    % --- Figure setup ---
    fig = figure('Name', 'Normalized Step Responses for Transfer Identification');
    fig.Position = [100 100 900 700];

    sgtitle({'Normalized Step Responses for System Identification', ...
        '$(y_i - y_{i,s}) / (u_j - u_{j,s})$', ...
        'Each subplot represents one transfer $G_{ij}(s)$'}, ...
        'Interpreter','latex', 'FontSize', 13);

    % --- Color and axes ---
    colors = lines(4);
    ax = gobjects(2,2);

    % === G11: h1 ← F1 ===
    ax(1,1) = subplot(2,2,1); hold on; grid on;
    plot(t, y_norm_F1(1,:), 'Color', colors(1,:), 'LineWidth', 1.8);
    xline(198, 'r--');
    text(198, 0.5, '$T_1$', 'Interpreter','latex', ...
         'Rotation',90, 'VerticalAlignment','bottom', 'HorizontalAlignment','right', Color='red');
    yline(0.632*0.301, 'r--', Label='$0.632*K_1$', Interpreter='latex')
    xlabel('Time [s]');
    ylabel('Normalized Response [(cm)/(cm$^3$/s)]', 'Interpreter','latex');
    title('$G_{11}$ : $h_1 \leftarrow F_1$', 'Interpreter','latex', 'FontWeight','bold');
    xlim([t(1), t(end)]);
    ylim([0, 0.5])

    % === G12: h1 ← F2 ===
    ax(1,2) = subplot(2,2,2); hold on; grid on;
    plot(t, y_norm_F2(1,:), 'Color', colors(2,:), 'LineWidth', 1.8);
    xlabel('Time [s]');
    ylabel('Normalized Response [(cm)/(cm$^3$/s)]', 'Interpreter','latex');
    title('$G_{12}$ : $h_1 \leftarrow F_2$', 'Interpreter','latex', 'FontWeight','bold');
    xlim([t(1), t(end)]);
    ylim([0, 0.5])

    % === G21: h2 ← F1 ===
    ax(2,1) = subplot(2,2,3); hold on; grid on;
    plot(t, y_norm_F1(2,:), 'Color', colors(3,:), 'LineWidth', 1.8);
    xlabel('Time [s]');
    ylabel('Normalized Response [(cm)/(cm$^3$/s)]', 'Interpreter','latex');
    title('$G_{21}$ : $h_2 \leftarrow F_1$', 'Interpreter','latex', 'FontWeight','bold');
    xlim([t(1), t(end)]);
    ylim([0, 0.5])

    % === G22: h2 ← F2 ===
    ax(2,2) = subplot(2,2,4); hold on; grid on;
    plot(t, y_norm_F2(2,:), 'Color', colors(4,:), 'LineWidth', 1.8);
    xline(220, 'r--');
    text(220, 0.5, '$T_2$', 'Interpreter','latex', ...
         'Rotation',90, 'VerticalAlignment','bottom', 'HorizontalAlignment','right', Color='red');
    yline(0.632*0.415, 'r--', Label='$0.632*K_{22}$', Interpreter='latex')
    xlabel('Time [s]');
    ylabel('Normalized Response [(cm)/(cm$^3$/s)]', 'Interpreter','latex');
    title('$G_{22}$ : $h_2 \leftarrow F_2$', 'Interpreter','latex', 'FontWeight','bold');
    xlim([t(1), t(end)]);
    ylim([0, 0.5])

    % --- Style adjustments ---
    set(findall(fig,'-property','FontSize'),'FontSize',12);
    set(gcf, 'PaperPositionMode','auto');
end
