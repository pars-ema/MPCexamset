
%   Add path to constants and functions' definitions
run('tank_module.m')

%%  2.4 - Simulation

X_det = zeros(4,N);
X_stoch = zeros(4,N);
X_sde = zeros(4,N);

u_all = zeros(2,N);

z_det = zeros(2,N);
z_stoch = zeros(2,N);
z_sde = zeros(2,N);


%   Objective height values (cm)
z_obj = [1;1];
d_state = [0;0];

for k=1:N
    %   Obtain piecewise constant u updated by z values
    u_all(:,k) = input_flow_u(k);

    %   d noise generation: F3, F4 (normalized to up to 20 cm3/s)
    d = rand(2)*norm_d;
    d = d(1,:);

    % Calculate Derivative (using the deterministic f)
    dX_det = det_nonlinear_f(X_det(:, k), u_all(:,k), p);
    dX_stoch = stoch_nonlinear_f(X_stoch(:, k), u_all(:,k), d, p);
    [dX_sde,d_diff,dw_it] = sde_nonlinear_f(X_sde(:, k),u_all(:,k),p,d_state,d_mean,dt,dw,dw_it);

    d_state = d_state + d_diff;

    % Forward Euler Step: X(k+1) = X(k) + dX/dt * dt
    if k < N
        X_det(:, k+1) = X_det(:, k) + dX_det * dt;
        X_stoch(:, k+1) = X_stoch(:, k) + dX_stoch * dt;
        X_sde(:, k+1) = X_sde(:, k) + dX_sde * dt;
    end
    
    % Calculate Output (z)
    z_det(:, k) = det_nonlinear_h(X_det(:, k),p);
    z_stoch(:, k) = stoch_nonlinear_h(X_stoch(:, k),p,Rd);
    z_sde(:, k) = sde_nonlinear_h(X_sde(:, k),p,Rd);
end


%   Plotting results
figure();
subplot(2,1,1);
hold on;
plot(t, z_det(1,:), 'b-', 'LineWidth', 1.5);
h_stoch = plot(t, z_stoch(1,:), 'r--', 'LineWidth', 1);
h_sde = plot(t, z_sde(1,:), 'g:', 'LineWidth', 1);

alpha = 0.3; 
set(h_stoch, 'Color', [1 0 0 alpha]);
set(h_sde, 'Color', [0 1 0 alpha]);

yline(z_obj(1), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Objective h1');
title('Tank 1 Height Response: $h_1$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Height (cm)');
legend('Deterministic', 'Stochastic (Process Noise $d$)', 'SDE (Process + Sigma Noise)', 'Location', 'best');
grid on;
hold off;

subplot(2,1,2);
hold on;
plot(t, z_det(2,:), 'b-', 'LineWidth', 1.5);
h_stoch = plot(t, z_stoch(2,:), 'r--', 'LineWidth', 1);
h_sde = plot(t, z_sde(2,:), 'g:', 'LineWidth', 1);

alpha = 0.3; 
set(h_stoch, 'Color', [1 0 0 alpha]);
set(h_sde, 'Color', [0 1 0 alpha]);

yline(z_obj(2), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Objective h2');
title('Tank 2 Height Response: $h_2$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Height (cm)');
legend('Deterministic', 'Stochastic (Process Noise $d$)', 'SDE (Process + Sigma Noise)', 'Location', 'best');
grid on;
hold off;

figure();
hold on;
plot(t, u_all(1,:), 'DisplayName', 'Input u1 / u2');
title("Tanks' Inputs (F1,F2)", 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Flow Rate ([cm3/s])');
grid on;
hold off


function u = input_flow_u(t_curr)
    if t_curr < 200
        u = [500.0; 500.0];     % [cm3/s] for t in [0, 200)
    elseif t_curr < 500
        u = [0.0; 0.0];         % [cm3/s] for t in [200, 400)
    elseif t_curr < 2000
        u = [500.0; 500.0];
    else
        u = [20.0; 20.0];
    end
end