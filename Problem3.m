% --------------------------------------------------------------
% Parameters
% --------------------------------------------------------------
a1 = 1.2272;        %[cm2] Area of outlet pipe 1
a2 = 1.2272;        %[cm2] Area of outlet pipe 2
a3 = 1.2272;        %[cm2] Area of outlet pipe 3
a4 = 1.2272;        %[cm2] Area of outlet pipe 4
A1 = 380.1327;      %[cm2] Cross sectional area of tank 1
A2 = 380.1327;      %[cm2] Cross sectional area of tank 2
A3 = 380.1327;      %[cm2] Cross sectional area of tank 3
A4 = 380.1327;      %[cm2] Cross sectional area of tank 4
gamma1 = 0.45;      % Flow distribution constant. Valve 1
gamma2 = 0.40;      % Flow distribution constant. Valve 2
g = 981;            % [cm/s2] The acceleration of gravity
rho = 1.00;         % [g/cm3] Density of water
p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

a = [a1; a2; a3; a4];
A = [A1; A2; A3; A4];
gamma = [gamma1; gamma2];
% --------------------------------------------------------------

% Simulation scenario

t0 = 0.0;                       % [s] Initial time
t_final = 20*60;                % [s] Final time
m10 = 0.0;                      % [g] Liquid mass in tank 1 at time t0
m20 = 0.0;                      % [g] Liquid mass in tank 2 at time t0
m30 = 0.0;                      % [g] Liquid mass in tank 3 at time t0
m40 = 0.0;                      % [g] Liquid mass in tank 4 at time t0
F1 = 0;                         % [cm3/s] Flow rate from pump 1
F2 = 0;                         % [cm3/s] Flow rate from pump 2
x0 = [m10; m20; m30; m40];
u = [F1; F2];                   % Initially set to 0

% --------------------------------------------------------------


%%  2.4 - Simulation

X_det = zeros(4,N);
X_stoch = zeros(4,N);
X_sde = zeros(4,N);

u_all = zeros(2,N);
du_all = zeros(2,N);

z_det = zeros(2,N);
z_stoch = zeros(2,N);
z_sde = zeros(2,N);

%   Objective height values
z_obj = [50;50];
d_diff_ = [0;0];

for k=1:N
    %   Obtain piecewise constant u updated by z values
    if k == 1
        error = [0;0];
        
        u_all(:,k) = P_controller(error);
        %u_all(:,k) = PI_controller(error,[0;0]);

        %u_all(:,k) = PID_controller(error,[0;0],[0;0]);
        %du_all(:,k) = u_all(:,k);
    else
        %   The deterministic result was chosen for error computation
        error = z_obj - z_det(:,k-1);
        
        u_all(:,k) = P_controller(error);
        %u_all(:,k) = PI_controller(error,u_all(:,k-1));

        %u_all(:,k) = PID_controller(error,u_all(:,k-1),du_all(:,k-1));
        %du_all(:,k) = u_all(:,k) - u_all(:,k-1);
    end

    %   d noise generation: F3, F4 (normalized to up to 50 cm3/s)
    d = rand(2)*norm_d;
    d = d(1,:);

    % Calculate Derivative (using the deterministic f)
    dX_det = det_nonlinear_f(X_det(:, k), u_all(:,k), p);
    dX_stoch = stoch_nonlinear_f(X_stoch(:, k), u_all(:,k), d, p);
    [dX_sde,d_diff,dw_it] = sde_nonlinear_f(X_sde(:, k), u_all(:,k), p, d_diff_, d_mean,dt,dw,dw_it);

    d_diff_ = d_diff;

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
figure(1);
subplot(2,1,1);
hold on;
plot(t, z_det(1,:), 'b-', 'LineWidth', 1.5);
plot(t, z_stoch(1,:), 'r--', 'LineWidth', 1);
plot(t, z_sde(1,:), 'g:', 'LineWidth', 1);
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
plot(t, z_stoch(2,:), 'r--', 'LineWidth', 1);
plot(t, z_sde(2,:), 'g:', 'LineWidth', 1);
yline(z_obj(2), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Objective h2');
title('Tank 2 Height Response: $h_2$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Height (cm)');
legend('Deterministic', 'Stochastic (Process Noise $d$)', 'SDE (Process + Sigma Noise)', 'Location', 'best');
grid on;
hold off;

%{
figure(2);
hold on;
plot(t, u_all(1,:), 'DisplayName', 'Input u1 / u2');
title("Tanks' Inputs (F1,F2)", 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Flow Rate ([cm3/s])');
grid on;
hold off
%}

%%  Controllers' Logic: Implementation - Uncomment & comment for test
%   P,PI,PID Controllers available for use.

function u = P_controller(error)
    u = 5000 .* error;
end

function u = PI_controller(error,prev_u)
    u = 0.01 .* [prev_u(1) 0; 0 prev_u(2)] * error + 5000 .* error;
end

function u = PID_controller(error,prev_u,du)
    u = 0.01 .* [prev_u(1) 0; 0 prev_u(2)] * error + 5000 .* error + ...
        0.03 .* [du(1) 0; 0 du(2)] * error;
end