function [fig, ax] = plot_heights_2(y_meas, y_truth, T, p, r)
%PLOT_HEIGHTS_2 Plot measured vs true tank heights over time
%   [fig, ax] = plot_heights_2(y_meas, y_truth, T, p, r)
%   y_meas : measured heights (4×N)
%   y_truth: true heights (4×N)
%   T      : time vector
%   p      : parameters (unused here but kept for consistency)
%   r      : (optional) reference for tanks 1 and 2

if nargin < 5
    r = [];  % no reference provided
end

order = [3, 4, 1, 2];
h_max = max(y_meas, [], "all");

% --- Create figure and axes handles ---
fig = figure;
sgtitle('Measured heights with noise');
ax = gobjects(1, 4); % Preallocate axes handles

for i = 1:4
    ax(i) = subplot(2, 2, order(i));
    plot(T, y_meas(i,:), 'LineWidth', 1.5, 'DisplayName', 'Measured'); hold on;
    plot(T, y_truth(i,:), 'r--', 'LineWidth', 2, 'DisplayName', 'True');

    % Plot reference for tanks 1 and 2 if provided
    if ~isempty(r)
        if i <= size(r,1)
            plot(T, r(i,:), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
        end
    end

    grid on;
    xlabel('Time [s]');
    ylabel(['h_' num2str(i) ' [cm]']);
    title(['Tank ' num2str(i)]);
    legend('show', 'Location', 'best');
    xlim([T(1), T(end)]);
    ylim([0, 1.1*h_max]);
end

end
