function [fig, ax] = plot_masses(X, T, p)
% Plot masses in the four tanks and return figure + axes handles

order = [3, 4, 1, 2];

% Convert to kg
X = X / 1000;
m_max = max(X, [], "all");

fig = figure; % create figure and store handle
sgtitle("Masses in the tanks");

ax = gobjects(1,4); % preallocate axes handles

for i = 1:4
    ax(i) = subplot(2,2,order(i)); % save each subplot handle
    plot(T, X(i,:), 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel(['m_' num2str(i) ' [kg]']);
    title(['Tank ' num2str(i)]);
    xlim([T(1) T(end)]);
    ylim([0 1.1*m_max]);
end

end
