
%   Add path to constants and functions' definitions
run('tank_module.m')

% Steady-State Calculation (for Open-Loop Step Experiments)
% --------------------------------------------------------------

F_ss = 20;
u_ss = [F_ss; F_ss];

% Calculate steady-state heights (h_ss = [h1; h2; h3; h4])
h3_ss = (1/(2*g)) * (((1-gamma2) * F_ss) / a3)^2;
h4_ss = (1/(2*g)) * (((1-gamma1) * F_ss) / a4)^2;
Q1_in_ss = gamma1*F_ss + a3*sqrt(2*g*h3_ss);
Q2_in_ss = gamma2*F_ss + a4*sqrt(2*g*h4_ss);
h1_ss = (1/(2*g)) * (Q1_in_ss / a1)^2;
h2_ss = (1/(2*g)) * (Q2_in_ss / a2)^2;
h_ss = [h1_ss; h2_ss; h3_ss; h4_ss];

X_ss = rho .* A .* h_ss;        % X_ss is the steady-state mass vector
z_ss = h_ss(1:2);               % steady-state output heights (h1, h2)

%   ----------------------------------------------------------------

%%  Continuous-time model definition: Deterministic Nonlinear Model

A_cd =      [   -((p(12)*p(1))/(2*sqrt(X_ss(1))))*sqrt(2*p(11)/(p(12)*p(5))) 0 ((p(12)*p(3))/(2*sqrt(X_ss(3))))*sqrt(2*p(11)/(p(12)*p(7))) 0;
                0 -((p(12)*p(2))/(2*sqrt(X_ss(2))))*sqrt(2*p(11)/(p(12)*p(6))) 0 ((p(12)*p(4))/(2*sqrt(X_ss(4))))*sqrt(2*p(11)/(p(12)*p(8)));
                0 0 -((p(12)*p(3))/(2*sqrt(X_ss(3))))*sqrt(2*p(11)/(p(12)*p(7))) 0;
                0 0 0 -((p(12)*p(4))/(2*sqrt(X_ss(4))))*sqrt(2*p(11)/(p(12)*p(8)))];

B_cd =      [   p(12)*p(9) 0;
                0 p(12)*p(10);
                0 p(12)*(1-p(10));
                p(12)*(1-p(9)) 0];

C_cd =      [   1/(p(12)*p(5)) 0 0 0;
                0 1/(p(12)*p(6)) 0 0];


%%  Continuous-time model definition: Stochastic Nonlinear Model

%   (F3,F4) do not make a difference in any of the matrices below. As 'u',
%   the vector [F1;F2] was considered, and [F3;F4] would be the disturbance
%   'd' vector. As their steady value would be (d_max - 0)/2, then only the
%   E matrix would be affected by them.

A_cs = A_cd;

B_cs = B_cd;

E_cs = [0 0;
        0 0;
        p(12) 0;
        0 p(12)];

C_cs = C_cd;


%%  Continuous-time model definition: Stochastic Nonlinear Model SDE

A_csde = A_cd;
B_csde = B_cd;
E_csde = E_cs;
C_csde = C_cd;

%   -----------------------------------------------------------------

%%  Problem 5.1: Compute continuous-time linearized models for the 3 
%   models developed in Problem 2.

X_det = zeros(4,N);
X_stoch = zeros(4,N);
X_sde = zeros(4,N);

u_all = zeros(2,N);

z_det = zeros(2,N);
z_stoch = zeros(2,N);
z_sde = zeros(2,N);

%   Objective height values (cm)
d_state = [0;0];
d_diff = [0;0];

for k=1:N
    %   Obtain piecewise constant u updated by z values
    u_all(:,k) = input_flow_u(k);

    %   d noise generation: F3, F4 (normalized to up to 20 cm3/s)
    d = rand(2)*norm_d;
    d = d(:,1);

    %   Calculate linear derivative
    dX_det = A_cd*(X_det(:, k) - X_ss) + B_cd*(u_all(:,k) - u_ss);
    dX_stoch = A_cs*(X_det(:, k) - X_ss) + B_cs*(u_all(:,k) - u_ss) + ...
                E_cs*(d - d_mean);
    dX_sde = A_csde*(X_det(:, k) - X_ss) + B_csde*(u_all(:,k) - u_ss) + ...
                E_csde*(d_state - d_mean);

    %   Update d_state for next SDE iteration
    [d_diff,~] = stoch_update_sigma(d_state,d_mean,dt,dw,dw_it);
    d_state = d_state + d_diff;

    %   Forward Euler Step: X(k+1) = X(k) + dX/dt * dt
    if k < N
        X_det(:, k+1) = X_det(:, k) + dX_det * dt;
        X_stoch(:, k+1) = X_stoch(:, k) + dX_stoch * dt;
        X_sde(:, k+1) = X_sde(:, k) + dX_sde * dt;
    end
    
    %   Calculate linear output (z)
    z_det(:, k) = C_cd*(X_det(:, k) - X_ss) + z_ss;
    z_stoch(:, k) = C_cs*(X_stoch(:, k) - X_ss) + z_ss + [normrnd(0, Rd(1));normrnd(0, Rd(2))];
    z_sde(:, k) = C_csde*(X_sde(:, k) - X_ss) + z_ss + [normrnd(0, Rd(1));normrnd(0, Rd(2))];
end


plot_results(z_det,z_stoch,z_sde,u_all,t)


%%  Problem 5.2: Compute the continuous-time transfer functions for the
%   continuous-time linearized models

%   Define the D matrix (feedthrough matrix, zero for this system)
D_cd = zeros(size(C_cd, 1), size(B_cd, 2));

%   Obtain G(s)
sys_ss = ss(A_cd, B_cd, C_cd, D_cd);
G_s = tf(sys_ss);                       %   G(s) matrix obtained

% 1. Calculate the input perturbation: Delta_u = u - u_ss
delta_u = u_all - u_ss; 

% 2. Transpose the input for lsim: 
% lsim requires input data to be [Time_steps x Inputs], i.e., [N x 2]
u_input = delta_u';

% Run the simulation: The output is the perturbation delta_z
[delta_z_sim, t_out] = lsim(G_s, u_input, t);

% Transpose back to (2xN)
delta_z_simt = delta_z_sim'; 

% Add the steady-state output back: z = Delta_z + z_ss
% z_sim is now a 2xN matrix of simulated heights (h1 and h2)
z_sim = delta_z_simt + z_ss;

%   Plotting results
figure();
subplot(2,1,1);
hold on;
plot(t, z_det(1,:), 'b-', 'LineWidth', 1.5);

title('Tank 1 Height Response G(s): $h_1$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Height (cm)');
legend('Deterministic G(s)', 'Location', 'best');
grid on;
hold off;

subplot(2,1,2);
hold on;
plot(t, z_det(2,:), 'b-', 'LineWidth', 1.5);

title('Tank 2 Height Response G(s): $h_2$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Height (cm)');
legend('Deterministic G(s)', 'Location', 'best');
grid on;
hold off;

K = -C_cd / A_cd * B_cd;
disp('Steady state gains:');
disp(K);

lambda = eig(A_cd);
tau = -1 ./ real(lambda);
disp('Time constants:');
disp(tau);


%%  Problem 5.3: Compare the gains and time constants to the gains and 
%   time constants obtained from the step response experiments in Problem 4.

%{
The tau constants are a little smaller than in the nonlinear approach given
from problem 4. But not so further away. And surprisingly, the K gains re-
sult quite similar.
%}


%%  Problem 5.4: Compute discrete-time state space models using a sampling
%   time of your choice (in this case, T_s = dt)

%   Obtain discrete matrices for ZOH inputs
sys_d = c2d(sys_ss,dt,'zoh');

%   Extract matrices
A_d = sys_d.A;
B_d = sys_d.B;
C_d = sys_d.C;
D_d = sys_d.D;

% Initialize state and output vectors for the discrete model
X_d = zeros(4, N);
z_d = zeros(2, N);

%   Discrete-Time Simulation using A_d, B_d
for k = 1:N
    %   Get delta_u
    delta_u = u_all(:, k) - u_ss;

    delta_z = C_d * (X_d(:, k) - X_ss) + D_d * delta_u;
    z_d(:, k) = delta_z + z_ss;

    %   Calculate next state (X[k+1]) from the current state (X[k])
    if k < N
        % X[k+1] - X_ss = A_d * (X[k] - X_ss) + B_d * U[k]
        X_d(:, k+1) = X_ss + A_d * (X_d(:, k) - X_ss) + B_d * delta_u;
    end
end

%   Plotting results
figure();
subplot(2,1,1);
hold on;
plot(t, z_d(1,:), 'b-', 'LineWidth', 1.5);

title('Tank 1 Height Response (discrete): $h_1$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Height (cm)');
legend('Linear (discrete)', 'Location', 'best');
grid on;
hold off;

subplot(2,1,2);
hold on;
plot(t, z_d(2,:), 'b-', 'LineWidth', 1.5);

title('Tank 2 Height Response (discrete): $h_2$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Height (cm)');
legend('Linear (discrete)', 'Location', 'best');
grid on;
hold off;


%%  Problem 5.5: Calculate the Markov Parameters

H = zeros(size(C_d,1),size(B_d,2),N);
factor = eye(size(A_d));

for k = 1:N
    if k == 1
        H(:,:,k) = D_d;
    else
        H(:,:,k) = C_d * factor * B_d;
        factor = factor * A_d;
    end
end


%   Plotting results
figure();
subplot(2,2,1);
hold on;
plot(t, squeeze(H(1,1,:)), 'LineWidth', 1.5);

title('Markov Parameter $H_1$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Magnitude ($s/cm^2$)','Interpreter', 'latex');
grid on;
hold off;

subplot(2,2,2);
hold on;
plot(t, squeeze(H(1,2,:)), 'LineWidth', 1.5);

title('Markov Parameter $H_2$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Magnitude ($s/cm^2$)','Interpreter', 'latex');
grid on;
hold off;

subplot(2,2,3);
hold on;
plot(t, squeeze(H(2,1,:)), 'LineWidth', 1.5);

title('Markov Parameter $H_3$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Magnitude ($s/cm^2$)','Interpreter', 'latex');
grid on;
hold off;

subplot(2,2,4);
hold on;
plot(t, squeeze(H(2,2,:)), 'LineWidth', 1.5);

title('Markov Parameter $H_4$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Magnitude ($s/cm^2$)','Interpreter', 'latex');
grid on;
hold off;


%%  Discuss and comment on the linearization approach for obtaining
%   discrete-time linear state space models.

%{
A pretty straightforward approach which brings simplicity for the system
and for the controlling process, along with pretty acceptable results if
an much more absolute precision is not required. 


%}


%%  Input flow u function, and plot function

function u = input_flow_u(t_curr)
    if t_curr < 200
        u = [500.0; 500.0];     % [cm3/s] for t in [0, 200)
    elseif t_curr < 500
        u = [0.0; 0.0];         % [cm3/s] for t in [200, 1000)
    elseif t_curr < 2000
        u = [500.0; 500.0];
    else
        u = [20.0; 20.0];       %   Steady state input
    end
end


function plot_results(z_det,z_stoch,z_sde,u_all,t)
    %   Plotting results
    figure();
    subplot(2,1,1);
    hold on;
    plot(t, z_det(1,:), 'b-', 'LineWidth', 1.5);
    h_stoch = plot(t, z_stoch(1,:), 'r--', 'LineWidth', 1);
    h_sde = plot(t, z_sde(1,:), 'g:', 'LineWidth', 1);
    
    alpha = 0.2; 
    set(h_stoch, 'Color', [1 0 0 alpha]);
    set(h_sde, 'Color', [0 1 0 alpha]);
    
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
    
    alpha = 0.2; 
    set(h_stoch, 'Color', [1 0 0 alpha]);
    set(h_sde, 'Color', [0 1 0 alpha]);
    
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
end
