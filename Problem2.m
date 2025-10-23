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
t_final = 10*60;                % [s] Final time
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

dt = 0.001;
t = t0:dt:t_final;
N = length(t);
norm_d = 15;

X_det = zeros(4,N);
X_stoch = zeros(4,N);
X_sde = zeros(4,N);

u_all = zeros(2,N);

z_det = zeros(2,N);
z_stoch = zeros(2,N);
z_sde = zeros(2,N);

Rd = [0.2 0.1 0.1 0.2];

%   Scalar wiener noise generation and iteration index
[~,~,dw] = ScalarStdWienerProcess(t(end),N,2,5);
dw_it = 1;

d_diff = [0.0; 0.0];
d_mean = [15.0;15.0];       %   Set to the max value of d

%   Objective height values (cm)
z_obj = [1;1];

for k=1:N
    %   Obtain piecewise constant u updated by z values
    u_all(:,k) = input_flow_u(k);

    %   d noise generation: F3, F4 (normalized to up to 20 cm3/s)
    d = rand(2)*norm_d;
    d = d(1,:);

    % Calculate Derivative (using the deterministic f)
    dX_det = det_nonlinear_f(X_det(:, k), u_all(:,k), p);
    dX_stoch = stoch_nonlinear_f(X_stoch(:, k), u_all(:,k), d, p);
    [dX_sde,d_diff,dw_it] = sde_nonlinear_f(X_sde(:, k),u_all(:,k),d,p,d_diff,d_mean,dt,dw,dw_it);

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

figure(2);
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
    elseif t_curr < 400
        u = [0.0; 0.0];         % [cm3/s] for t in [200, 400)
    elseif t_curr < 2000
        u = [500.0; 500.0];
    else
        u = [60.0; 60.0];
    end
end


%%  2.1 - Deterministic Nonlinear Model development

%   Prediction function f(X,u,p)
%   x = [m1; m2; m3; m4] (masses in g)
%   u = [F1; F2] (pump flow rates in cm3/s)
function dX = det_nonlinear_f(X,u,p)
    dX = [0;0;0;0];
    dX(1) =  p(12)*(p(9)*u(1)) + p(12)*p(3)*sqrt(2*p(11)*(X(3)/(p(12)*p(7)))) ...
            - p(12)*p(1)*sqrt(2*p(11)*(X(1)/(p(12)*p(5))));
    dX(2) =  p(12)*(p(10)*u(2)) + p(12)*p(4)*sqrt(2*p(11)*(X(4)/(p(12)*p(8)))) ...
            - p(12)*p(2)*sqrt(2*p(11)*(X(2)/(p(12)*p(6))));
    dX(3) =  p(12)*((1-p(10))*u(2)) - p(12)*p(3)*sqrt(2*p(11)*(X(3)/(p(12)*p(7))));
    dX(4) =  p(12)*((1-p(9))*u(1)) - p(12)*p(4)*sqrt(2*p(11)*(X(4)/(p(12)*p(8))));
end

%   Measurement function g(X,p)
function h = det_nonlinear_g(X,p)
    % Height h_i = m_i / (rho * A_i)
    h = [0;0;0;0];
    h(1) = X(1) / (p(12) * p(5));
    h(2) = X(2) / (p(12) * p(6));
    h(3) = X(3) / (p(12) * p(7));
    h(4) = X(4) / (p(12) * p(8));
    
end

%   The heights (h_x) are the outputs
function z = det_nonlinear_h(X,p)
    z = det_nonlinear_g(X,p);
    z = z(1:2);
end



%%  2.2 - Stochastic Nonlinear Model development

%   (F3,F4) are now defined and included in the scope of this exercise.
%   Those should be defined stochastic but piecewise constant

%   d = [F3; F4] (stochastic flow disturbances in cm3/s)
function dX = stoch_nonlinear_f(X,u,d,p)
    dX = det_nonlinear_f(X,u,p);
    
    dX(3) = dX(3) + d(1);
    dX(4) = dX(4) + d(2);
end


%   2. Measurement function g(x,p)
%   The measurements are the liquid heights (h_i) in cm.
%   y(t) = g(x(t),p) + v(t)
function h = stoch_nonlinear_g(X,p,R)
    %   For this application, R should be a vector 1x4, representing a
    %   diagonal matrix.
    v = rand(length(R));
    v = R*v;

    % Height h_i = m_i / (rho * A_i)
    h = [0;0;0;0];
    h(1) = (X(1) / (p(12) * p(5))) + v(1);
    h(2) = (X(2) / (p(12) * p(6))) + v(2);
    h(3) = (X(3) / (p(12) * p(7))) + v(3);
    h(4) = (X(4) / (p(12) * p(8))) + v(4);
    
end

function z = stoch_nonlinear_h(X,p,R)
    z = stoch_nonlinear_g(X,p,R);
    z = z(1:2);
end



%%  2.3 - Stochastic Nonlinear Model SDE


function [W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed)
    % ScalarStdWienerProcess Ns realizations of a scalar std Wiener process
    %
    % Syntax: [W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed)
    % W : Standard Wiener process in [0,T]
    % Tw : Time points
    % dW : White noise used to generate the Wiener process
    %
    % T : Final time
    % N : Number of intervals
    % Ns : Number of realizations
    % seed : To set the random number generator (optional)
    if nargin == 4
        rng(seed);
    end

    dt = T/N;
    dW = sqrt(dt)*randn(Ns,N);
    W = [zeros(Ns,1) cumsum(dW,2)];
    Tw = 0:dt:T;
end


%   Derivation for the stochastic added variables F3 and F4
function [d_diff,dw_it] = stoch_update_sigma(d,d_mean,dt,dw,dw_it)
    tau_v = [0.2; 0.1];
    sigma_d = [0.1; 0.1];

    d_diff = tau_v.*(d_mean - d)*dt + sigma_d.*(dw(:,dw_it));
    dw_it = dw_it + 1;
end

%   The sigma noise is given as a diagonal matrix/format
function [dX,d_diff,dw_it] = sde_nonlinear_f(X,u,d,p,d_diff,d_mean,dt,dw,dw_it)
    dX = stoch_nonlinear_f(X,u,d_diff,p);
    [d_diff,dw_it] = stoch_update_sigma(d,d_mean,dt,dw,dw_it);
end

function h = sde_nonlinear_g(X,p,R)
    h = stoch_nonlinear_g(X,p,R);
end

function z = sde_nonlinear_h(X,p,R)
    z = sde_nonlinear_g(X,p,R);
    z = z(1:2);
end