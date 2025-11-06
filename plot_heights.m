function [fig, ax] = plot_heights(X, T, p, r)
%PLOT_HEIGHTS Plot tank heights over time
%   [fig, ax] = plot_heights(X, T, p, r)
%   X = states (mass)
%   T = time vector
%   p = parameters
%   r = (optional) reference for tanks 1 and 2

if nargin < 4
    r = [];  % no reference provided
end

order = [3, 4, 1, 2];

A = p(5:8,1);   % Tank cross-sectional areas
rho = p(12);    % Density

[nX, nT] = size(X);
H = zeros(nX, nT);

% Compute heights from masses
for i = 1:nT
    H(:,i) = X(:,i) ./ (rho * A);
end

h_max = max(H, [], "all");

% --- Create figure and axes ---
fig = figure;
sgtitle('Real heights based on masses');
ax = gobjects(1, 4); % Preallocate axes handles

for i = 1:4
    ax(i) = subplot(2,2,order(i));
    plot(T, H(i,:), 'LineWidth', 1.5); hold on;

    % Plot reference for tanks 1 and 2 if provided
    if ~isempty(r)
        if i == 1
            plot(T, r(1)*ones(size(T)), 'k--', 'LineWidth', 1.5);
        elseif i == 2
            plot(T, r(2)*ones(size(T)), 'k--', 'LineWidth', 1.5);
        end
    end

    grid on;
    xlabel('Time [s]');
    ylabel(['h_' num2str(i) ' [cm]']);
    title(['Tank ' num2str(i)]);
    xlim([T(1) T(end)]);
    ylim([0 1.1*h_max]);
end

end
