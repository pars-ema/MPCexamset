function graph = plot_noisy_meassurements(X, T, p, show_truth, r)
%PLOT_HEIGHTS Plot tank heights over time
%   X = states (mass)
%   T = time vector
%   p = parameters
%   r = (optional) reference for tanks 1 and 2

if nargin < 5
    r = [];  % no reference provided
end

order = [3, 4, 1, 2];

[nX,nT] = size(X);
H = zeros(nX,nT);

if show_truth == 1
    H_truth = zeros(nX,nT);
end

N = length(T);

R = eye(4);
Lr = chol(R,"lower");
v = Lr*randn(4,N)



% Compute heights from masses
for i = 1:nT
    H(:, i) = FourTankSystemSensor(X(:, i) ,p) + v(:, i);
    if show_truth == 1
        H_truth(:, i) = FourTankSystemSensor(X(:, i) ,p);
    end
end



h_max = max(H,[], "all");

figure;
sgtitle('Noisy meassurements of heights'); 

for i = 1:4
    subplot(2,2,order(i));
    plot(T, H(i,:), 'LineWidth', 1.5); hold on;

    if show_truth == 1
        plot(T, H_truth(i,:), 'LineWidth', 1.5, LineStyle='--', Color='red')
    end



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

graph = 1;
end