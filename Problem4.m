% Four-Tank System: Step Response Simulation and Identification

% --------------------------------------------------------------
% Parameters (from original script)
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
A = [A1; A2; A3; A4];

% --------------------------------------------------------------
% Steady-State Calculation (for Open-Loop Step Experiments)
% --------------------------------------------------------------

F_ss = 50;
u_ss = [F_ss; F_ss];

% Calculate steady-state heights (h_ss = [h1; h2; h3; h4])
h3_ss = (1/(2*g)) * (((1-gamma2) * F_ss) / a3)^2;
h4_ss = (1/(2*g)) * (((1-gamma1) * F_ss) / a4)^2;
Q1_in_ss = gamma1*F_ss + a3*sqrt(2*g*h3_ss);
Q2_in_ss = gamma2*F_ss + a4*sqrt(2*g*h4_ss);
h1_ss = (1/(2*g)) * (Q1_in_ss / a1)^2;
h2_ss = (1/(2*g)) * (Q2_in_ss / a2)^2;
h_ss = [h1_ss; h2_ss; h3_ss; h4_ss];

x0 = rho .* A .* h_ss;      % x0 is the steady-state mass vector
z_ss = h_ss(1:2);           % steady-state output heights (h1, h2)

% --------------------------------------------------------------
% Simulation Parameters
% --------------------------------------------------------------
t_final = 1200;                             % [s] Final time (20 min)
dt = 1;                                     % [s] Sampling time
t = 0:dt:t_final;
N = length(t);
t_step = 60;                                % Step time at 1 minute
k_step = find(t >= t_step, 1);
step_percentages = [0.10, 0.25, 0.50];
step_mags = F_ss * step_percentages;        % Step magnitude (F_step - F_ss)
nu = 2;                                     % Number of inputs (F1, F2)
ny = 2;                                     % Number of outputs (h1, h2)

% --------------------------------------------------------------


%% 1. Simulate the step responses for 10%, 25% and 50% steps (Deterministic Model)
% Store results: (output_dim, time, step_size_index, input_tested_index)
Z_det = zeros(ny, N, length(step_percentages), nu); 
for i = 1:nu
    for j = 1:length(step_percentages)
        step_mag = step_mags(j);
        
        X_det = zeros(4, N);
        X_det(:, 1) = x0;
        u_all = repmat(u_ss, 1, N);
        
        % Apply step to the i-th input
        u_all(i, k_step:end) = u_ss(i) + step_mag;
       
        % Calculate Output z = [h1; h2] and X vector
        for k = 1:N
            if k < N
                dX_det = det_nonlinear_f(X_det(:, k), u_all(:,k), p);
                X_det(:, k+1) = X_det(:, k) + dX_det * dt;
            end
            
            Z_det(:, k, j, i) = det_nonlinear_h(X_det(:, k), p);
        end
    end
end


%% 2. Simulate the step responses with process and measurement noise

noise_levels = {'low', 'medium', 'high'};
noise_factors = [0.1, 0.3, 0.6];
base_Rd = [0.2, 0.1, 0.1, 0.2];
max_disturb = 5;

% Store results: (output_dim, time, step_size_index, input_tested_index, noise_level_index)
Z_stoch = zeros(ny, N, length(step_percentages), nu, length(noise_levels));
for k_noise = 1:length(noise_levels)
    Rd = base_Rd * noise_factors(k_noise);              % Measurement Noise Vector
    norm_d = max_disturb * noise_factors(k_noise);      % Max Disturbance for F3, F4
    for i = 1:nu
        for j = 1:length(step_percentages)
            step_mag = step_mags(j);
            
            X_stoch = zeros(4, N);
            X_stoch(:, 1) = x0;
            u_all = repmat(u_ss, 1, N);
            u_all(i, k_step:end) = u_ss(i) + step_mag;
            % Simulation loop
            for k = 1:N-1
                % Process Noise 'd' (F3, F4) - Uniformly distributed from [0, norm_d]
                d = rand(2, 1) * norm_d; 
                
                % Stochastic Process: stoch_nonlinear_f
                dX_stoch = stoch_nonlinear_f(X_stoch(:, k), u_all(:,k), d, p);
                X_stoch(:, k+1) = X_stoch(:, k) + dX_stoch * dt;
                
                % Measurement Noise: stoch_nonlinear_h
                Z_stoch(:, k, j, i, k_noise) = stoch_nonlinear_h(X_stoch(:, k), p, Rd);
            end
            % Final output (k=N) calculation
            Z_stoch(:, N, j, i, k_noise) = stoch_nonlinear_h(X_stoch(:, N), p, Rd);
        end
    end
end


%% 3. Compute and plot (in appropriate plots) the normalized steps

% 3a. Deterministic Normalized Steps
Z_norm_det = zeros(ny, N, length(step_percentages), nu);
figure(1);
sgtitle('3a. Normalized Step Responses (Deterministic)');
for i = 1:nu
    for l = 1:ny
        subplot(nu, ny, (i-1)*ny + l);
        hold on;
        for j = 1:length(step_percentages)
            step_mag = step_mags(j); 

            % Normalized Response = (h(t) - h_ss) / Delta_u
            Z_norm_det(l, :, j, i) = (Z_det(l, :, j, i) - z_ss(l)) / step_mag; 
            plot(t/60, Z_norm_det(l, :, j, i), 'LineWidth', 1.5, 'DisplayName', ...
                sprintf('%d%% Step', step_percentages(j)*100));
        end
        title(sprintf('Normalized Step: $\\Delta h_{%d}$ / $\\Delta F_{%d}$', l, i), ...
                'Interpreter', 'latex');
        xlabel('Time (min)');
        ylabel('Normalized Height (cm / (cm^3/s))');
        legend('Location', 'best');
        grid on;
        hold off;
    end
end

%   Z_norm_stoch now stores all results:
%   (output_dim, time, step_size_index, input_tested_index, noise_level_index)
Z_norm_stoch = zeros(ny, N, length(step_percentages), nu, length(noise_levels));

% Calculate ALL Normalized Stochastic Responses
for k_noise = 1:length(noise_levels)
    for i = 1:nu % Inputs (F1, F2)
        for j = 1:length(step_percentages)
            step_mag = step_mags(j);

            % Normalized Response = (h(t) - h_ss) / Delta_u
            Z_norm_stoch(:, :, j, i, k_noise) = (Z_stoch(:, :, j, i, k_noise) - z_ss) / step_mag;
        end
    end
end

% Plotting - Compare Step Sizes (Fixed Noise: Low)
k_noise_level = 1; % Low noise
figure(2);
sgtitle(sprintf('Stochastic Normalized Steps - Varying Step Size (%s Noise)', noise_levels{k_noise_level}));
for i = 1:nu
    for l = 1:ny % Outputs (h1, h2)
        subplot(nu, ny, (i-1)*ny + l);
        hold on;
        for j = 1:length(step_percentages)
            plot(t/60, Z_norm_stoch(l, :, j, i, k_noise_level), 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Stoch. %d%% Step', step_percentages(j)*100));
            plot(t/60, Z_norm_det(l, :, j, i), '--k', 'LineWidth', 0.5, ...
                 'HandleVisibility', 'off');
        end
        title(sprintf('$\\Delta h_{%d}$ / $\\Delta F_{%d}$', l, i), 'Interpreter', 'latex');
        xlabel('Time (min)');
        ylabel('Normalized Height');
        legend('Location', 'best');
        grid on;
        hold off;
    end
end


% Plotting - Compare Step Sizes (Fixed Noise: Medium)
k_noise_level = 2; % Medium noise
figure(3);
sgtitle(sprintf('Stochastic Normalized Steps - Varying Step Size (%s Noise)', noise_levels{k_noise_level}));
for i = 1:nu
    for l = 1:ny % Outputs (h1, h2)
        subplot(nu, ny, (i-1)*ny + l);
        hold on;
        for j = 1:length(step_percentages)
            plot(t/60, Z_norm_stoch(l, :, j, i, k_noise_level), 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Stoch. %d%% Step', step_percentages(j)*100));
            plot(t/60, Z_norm_det(l, :, j, i), '--k', 'LineWidth', 0.5, ...
                 'HandleVisibility', 'off');
        end
        title(sprintf('$\\Delta h_{%d}$ / $\\Delta F_{%d}$', l, i), 'Interpreter', 'latex');
        xlabel('Time (min)');
        ylabel('Normalized Height');
        legend('Location', 'best');
        grid on;
        hold off;
    end
end


% Plotting - Compare Step Sizes (Fixed Noise: High)
k_noise_level = 3; % High noise
figure(4);
sgtitle(sprintf('Stochastic Normalized Steps - Varying Step Size (%s Noise)', noise_levels{k_noise_level}));
for i = 1:nu
    for l = 1:ny % Outputs (h1, h2)
        subplot(nu, ny, (i-1)*ny + l);
        hold on;
        for j = 1:length(step_percentages)
            plot(t/60, Z_norm_stoch(l, :, j, i, k_noise_level), 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Stoch. %d%% Step', step_percentages(j)*100));
            plot(t/60, Z_norm_det(l, :, j, i), '--k', 'LineWidth', 0.5, ...
                 'HandleVisibility', 'off');
        end
        title(sprintf('$\\Delta h_{%d}$ / $\\Delta F_{%d}$', l, i), 'Interpreter', 'latex');
        xlabel('Time (min)');
        ylabel('Normalized Height');
        legend('Location', 'best');
        grid on;
        hold off;
    end
end

%{
Interpretation of the Plots (Nonlinearity)
The figure shows that the system is nonlinear because the
steady-state value (the normalized gain) is different for each input step size:

- **10% Step (Blue Line):** This line reaches the **highest** steady-state value.
    This is the best linear approximation, as the small step keeps the system
    closest to the initial operating point, where the linear model is most accurate.
- **25% Step (Orange Line):** The steady-state value is **lower** than the 10% step.
- **50% Step (Green Line):** The steady-state value is the **lowest**.

This is the result for the stochastic simulation of the system. Apparently
the K decreases when the delta(F) increases, but in fact that happens
because of the noise interference:

Response = (delta(h) + v) / delta(F)

For the deterministic plot, on the contrary, the 10% step magnitude signal
produces the lowest stationary output value. The gain increases as the
delta(F) increases. In this system there is no v.
Hence, in the previous system the decrease(F)-->increase(response)
potentially happens because of the influence of the v value, and the
relationship between delta(h) / delta(F).
%}


%% 4. From the normalized steps, identify a transfer function

%   Use the deterministic, 10% step normalized response for system 
%   identification, as this provides the cleanest, most linear approximation. 
%   We'll identify a First-Order Plus Dead Time (FOPDT) model for the 
%   2x2 MIMO system

% Re-initialize variables based on deterministic 10% step
k_step_id = 1;
nu = 2; ny = 2; N = 1201; dt = 1;

% Normalized response used for identification: Z_norm_id(output_i, time, input_j)
Z_norm_id = squeeze(Z_norm_det(:, :, k_step_id, :)); 

% 4a. Extract Final Gains ($K_{ij}$)
% K is the steady-state value of the normalized step response
% (average of last 10% of simulation).
K_ij = mean(Z_norm_id(:, round(0.9*N):N, :), 2);
K_ij = squeeze(K_ij);

K11_calc = K_ij(1, 1);      % ~0.0177
K12_calc = K_ij(1, 2);      % ~0.0244
K21_calc = K_ij(2, 1);      % ~0.0202
K22_calc = K_ij(2, 2);      % ~0.0142

% 4b. Estimated FOPDT Parameters (Tau_ij, Theta_ij)
% By inspection of the normalized step responses:
theta11 = 0;
theta12 = 0;
theta21 = 0;
theta22 = 0;

%   Calculate tau by 63.2% value of the stationary output

Target_11 = K11_calc * 0.632;
k_63_11 = find(Z_norm_id(1, :, 1) >= Target_11, 1, 'first') - 1;
t_63_11 = t(k_63_11);
tau11 = t_63_11 - t_step - theta11;

Target_22 = K22_calc * 0.632;
k_63_22 = find(Z_norm_id(2, :, 2) >= Target_22, 1, 'first') - 1;
t_63_22 = t(k_63_22);
tau22 = t_63_22 - t_step - theta22;

% --- Cross-Coupling Paths (G12, G21): Tangential Slope Method ---
% Method: tau = K / S, where S is the maximum slope (rate of change)
% The maximum slope is the maximum value of the numerical derivative (diff)
% of the response.

Response_12 = Z_norm_id(1, k_step:N, 2); 
Derivative_12 = diff(Response_12) / dt;
S_12 = max(abs(Derivative_12));
tau12 = abs(K12_calc) / S_12; 

Response_21 = Z_norm_id(2, k_step:N, 1);
Derivative_21 = diff(Response_21) / dt;
S_21 = max(abs(Derivative_21));
tau21 = abs(K21_calc) / S_21;

G11_s = sprintf('(%.3f e^{-%.0f s}) / (%.0f s + 1)', K11_calc, theta11, tau11);
G12_s = sprintf('(%.3f e^{-%.0f s}) / (%.0f s + 1)', K12_calc, theta12, tau12);
G21_s = sprintf('(%.3f e^{-%.0f s}) / (%.0f s + 1)', K21_calc, theta21, tau21);
G22_s = sprintf('(%.3f e^{-%.0f s}) / (%.0f s + 1)', K22_calc, theta22, tau22);

%   Create the Transfer Function Matrix G(s) using the identified FOPDT parameters
G11 = tf(K11_calc, [tau11, 1], 'InputDelay', theta11);
G12 = tf(K12_calc, [tau12, 1], 'InputDelay', theta12);
G21 = tf(K21_calc, [tau21, 1], 'InputDelay', theta21);
G22 = tf(K22_calc, [tau22, 1], 'InputDelay', theta22);

G_tf = [G11, G12; G21, G22]; % The 2x2 Transfer Function Matrix G(s)

disp(' ');
disp('4. Identified Transfer Function Matrix G(s) from u = [F1; F2] to y = [h1; h2]:');
disp('G(s) = [ G11(s)  G12(s) ]');
disp('       [ G21(s)  G22(s) ]');
disp(' ');
fprintf('[ %-50s %-50s ]\n', G11_s, G12_s);
fprintf('[ %-50s %-50s ]\n', G21_s, G22_s);


%% 5. Report the identified linear model estimate and discuss accuracy

%   Define Input for the linear model (Delta_U = 10% step)
k_step_id = 1;
Delta_u = step_mags(k_step_id);     %   Magnitude of the 10% step (F_ss * 0.1)

%   T1: Step applied only to F1 (U = [Delta_F1; 0])
U1 = zeros(nu, N);
U1(1, k_step:end) = Delta_u; 
Y_lin1 = lsim(G_tf, U1', t)';

%   T2: Step applied only to F2 (U = [0; Delta_F2])
U2 = zeros(nu, N);
U2(2, k_step:end) = Delta_u; 
Y_lin2 = lsim(G_tf, U2', t)';

%   Combine the results for plotting: Y_lin(output, time, input_tested_index)
Y_lin = zeros(ny, N, nu);
Y_lin(:, :, 1) = Y_lin1;
Y_lin(:, :, 2) = Y_lin2;

%   Calculate Nonlinear Output Change (h(t) - h_ss)
Y_nonlin = zeros(ny, N, nu);
for i = 1:nu
    Y_nonlin(:, :, i) = Z_det(:, :, k_step_id, i) - z_ss;
end

%   Plot the comparison
figure(5);
sgtitle('5. Model Accuracy: FOPDT Linear Model vs. Nonlinear System (50% Step)');

for i = 1:nu
    for l = 1:ny
        subplot(nu, ny, (i-1)*ny + l);
        hold on;
        
        plot(t/60, Y_nonlin(l, :, i), 'b-', 'LineWidth', 2, 'DisplayName', 'Nonlinear System (Actual)');
        plot(t/60, Y_lin(l, :, i), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Linear FOPDT Model (Prediction)');
        
        title(sprintf('$G_{%d,%d}$ Response ($\\Delta h_{%d}$ / $\\Delta F_{%d}$)', l, i, l, i), 'Interpreter', 'latex');
        xlabel('Time (min)');
        ylabel('Height Change $\\Delta h$ (cm)');
        legend('Location', 'best');
        grid on;
        hold off;
    end
end

%{
Discussion of Model Accuracy (Refer to Figure 5):
1. ACCURACY: The FOPDT model provides an excellent fit for the 10% step response.
   The small step ensures the system operates near the linearization point, minimizing
   errors from the nonlinear dynamics. There is some deviation on the
   cross-calculated delta(h_21)/delta(F_21) and delta(h_12)/delta(F_12),
   due to the lac of observability on the tanks 3 and 4.
2. NONLINEARITY: Accuracy will decrease significantly for larger steps (25%, 50%)
   because the true steady-state gain changes with input magnitude.
3. PHYSICAL INSIGHT: The model simulates correctly the first-order dynamics
   in the direct paths (G11, G22) and the cross-coupling paths (G12, G21)
   due to unmeasured tank delays.

Requirements of a Step Experiment:
1. INDEPENDENT STEPS: Inputs (F1 and F2) must be stepped **one at a time** to
   decouple the system and uniquely identify all four G(s) elements.
2. SMALL MAGNITUDE: The step must be **small (e.g., 10%)** to ensure the dynamics
   stay within the linear region for accurate parameter estimation.
3. DURATION: The test must run long enough (e.g., 1200s) to reach the new
   **steady-state** for accurate gain measurement.
%}


%% 6. Compute the corresponding impulse response coefficients (Markov Parameters)
% Sampling time is chosen as dt = 1s.

% Use the deterministic, 10% step normalized response (Z_norm_id) for M_k.
M = zeros(ny, N, nu); % M(output_i, time, input_j)

% Compute the difference (discrete-time derivative)
for i = 1:nu
    for l = 1:ny
        M(l, 2:N, i) = diff(Z_norm_id(l, :, i)) / dt;
    end
end

% Plotting the Markov Parameters
figure(6);
sgtitle('6. Impulse Response Coefficients (Markov Parameters)');
for i = 1:nu
    for l = 1:ny
        subplot(nu, ny, (i-1)*ny + l);
        hold on;

        % Plot starts from k=1 (i.e., time=dt)
        stem(t(2:N)/60, M(l, 2:N, i), 'filled', 'DisplayName', 'Markov Parameters'); 
        title(sprintf('$M_{%d,%d} (\\Delta h_{%d} / \\Delta F_{%d})$', l, i, l, i), 'Interpreter', 'latex');
        xlabel('Time (min)');
        ylabel('Coefficient Value');
        grid on;
        hold off;
    end
end


%{
The Markov Parameters (M's) start with a large positive value at t=0 and
then decay exponentially toward zero. It means that the system responds
immediately and very strongly to the input change (no dead time, θ=0). The
response quickly dissipates as the system settles. The initial large peak
represents the fastest charging rate of the tanks
%}


%%  Deterministic Nonlinear Model development

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



%%  Stochastic Nonlinear Model development

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



%%  Stochastic Nonlinear Model SDE

%   The sigma noise is given as a diagonal matrix/format
function dX = sde_nonlinear_f(X,u,d,p,sigma_diag)
    dX = stoch_nonlinear_f(X,u,d,p);

    dX(1) = dX(1) + sigma_diag(1);
    dX(2) = dX(2) + sigma_diag(2);
    dX(3) = dX(3) + sigma_diag(3);
    dX(4) = dX(4) + sigma_diag(4);
end

function h = sde_nonlinear_g(X,p,R)
    h = stoch_nonlinear_g(X,p,R);
end

function z = sde_nonlinear_h(X,p,R)
    z = sde_nonlinear_g(X,p,R);
    z = z(1:2);
end