function plot_disturbances(T, D)
% PLOT_DISTURBANCES Plot disturbance flows F3 and F4 over time
%
% Syntax: plot_disturbances(T, D)
%
% Inputs:
%   T - time vector [s]
%   D - matrix of disturbances [length(T) x 2]
%       D(:,1) = F3 (disturbance to tank 3)
%       D(:,2) = F4 (disturbance to tank 4)

figure;
stairs(T, D(1,:), 'b-', 'LineWidth', 1.5); hold on;
stairs(T, D(2,:), 'r-', 'LineWidth', 1.5);

xlabel('Time [s]');
ylabel('Disturbance flow [cm^3/s]');
title('Disturbance Signals F3 and F4');
legend('F3 (to tank 3)', 'F4 (to tank 4)');
grid on;
end
