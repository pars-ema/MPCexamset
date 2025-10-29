% --- Configuration and Setup ---
run('tank_module.m')

% --- Steady-State Calculation ---
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
x0 = rho .* A .* h_ss;      % Initial state (mass vector)
z_ss = h_ss(1:2);           % Steady-state output heights (h1, h2)

% --- Simulation Parameters ---
t_step = 60;                % Step time at 1 minute (60s)
% Find the time index where the step is applied.
k_step = find(t >= t_step, 1); 
step_percentages = [0.10, 0.25, 0.50];
nu = 2; ny = 2; 

noise_factors = [0.1 0.25 0.5];

% --- Process Noise Initialization (Fix for undefined 'norm_d') ---
max_disturb = 5;
norm_d_low = max_disturb * noise_factors(1);
norm_d_medium = max_disturb * noise_factors(2);
norm_d_high = max_disturb * noise_factors(3);
norm_d_list = [norm_d_low, norm_d_medium, norm_d_high];
step_mags = F_ss * step_percentages;

%% 1. Simulate the step responses (Deterministic Model)
Z_det = zeros(ny, N, length(step_percentages), nu); 
for i = 1:nu
    for j = 1:length(step_percentages)
        step_mag = step_mags(j);
        X_det = zeros(4, N);
        X_det(:, 1) = x0;
        u_all = repmat(u_ss, 1, N);
        u_all(i, k_step:end) = u_ss(i) + step_mag;
       
        for k = 1:N-1
            dX_det = det_nonlinear_f(X_det(:, k), u_all(:,k), p);
            X_det(:, k+1) = X_det(:, k) + dX_det * dt;
            
            % CRITICAL FIX: State Saturation (Mass/Height cannot be negative)
            X_det(:, k+1) = max(0, X_det(:, k+1)); 
            
            Z_det(:, k, j, i) = det_nonlinear_h(X_det(:, k), p);
        end
        Z_det(:, N, j, i) = det_nonlinear_h(X_det(:, N), p);
    end
end

%% 2. Simulate the step responses with noise (Stochastic Model)
Z_stoch = zeros(ny, N, length(step_percentages), nu, length(noise_factors));
for k_noise = 1:length(noise_factors)
    C_Rd = Rd * noise_factors(k_noise);
    norm_d = norm_d_list(k_noise); % Use the correctly defined norm_d
    
    for i = 1:nu
        for j = 1:length(step_percentages)
            step_mag = step_mags(j);
            X_stoch = zeros(4, N);
            X_stoch(:, 1) = x0;
            u_all = repmat(u_ss, 1, N);
            u_all(i, k_step:end) = u_ss(i) + step_mag;
            
            for k = 1:N-1
                % Process Noise 'd' (F3, F4) - Uniformly distributed from [0, norm_d]
                d = rand(2, 1) * norm_d; 
                
                dX_stoch = stoch_nonlinear_f(X_stoch(:, k), u_all(:,k),d,p);
                X_stoch(:, k+1) = X_stoch(:, k) + dX_stoch * dt;
                
                % CRITICAL FIX: State Saturation
                X_stoch(:, k+1) = max(0, X_stoch(:, k+1));
                
                % Measurement Noise: stoch_nonlinear_h
                Z_stoch(:, k, j, i, k_noise) = stoch_nonlinear_h(X_stoch(:, k), p, C_Rd);
            end
            Z_stoch(:, N, j, i, k_noise) = stoch_nonlinear_h(X_stoch(:, N), p, C_Rd);
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
Z_norm_stoch = zeros(ny, N, length(step_percentages), nu, length(noise_factors));

% Calculate ALL Normalized Stochastic Responses
for k_noise = 1:length(noise_factors)
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
sgtitle(sprintf('Stochastic Normalized Steps - Varying Step Size (%s Noise)', noise_factors(k_noise_level)));
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
sgtitle(sprintf('Stochastic Normalized Steps - Varying Step Size (%s Noise)', noise_factors(k_noise_level)));
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
sgtitle(sprintf('Stochastic Normalized Steps - Varying Step Size (%s Noise)', noise_factors(k_noise_level)));
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


%% 4. Automated Transfer Function Identification (Refactored & Simplified)
% Use the deterministic, 10% step normalized response (Z_norm_id) for system 
% identification, as this provides the cleanest, most linear approximation. 

k_step_id = 1;
Delta_u_id = step_mags(k_step_id);     % Magnitude of the 10% step

% Input Step for F1 (U = [Delta_F1; 0])
u_F1_step = repmat(u_ss, 1, N);
u_F1_step(1, k_step:end) = u_ss(1) + Delta_u_id; 

% The perturbation vector Delta_u_F1: Full signal minus steady state (u_ss)
Delta_u_F1 = u_F1_step - u_ss;         

% Output Perturbation Delta_z_F1: Output signal minus steady state (z_ss)
Z_F1 = Z_det(:, :, k_step_id, 1);
Delta_z_F1 = Z_F1 - z_ss;

% Input Step for F2 (U = [0; Delta_F2])
u_F2_step = repmat(u_ss, 1, N);
u_F2_step(2, k_step:end) = u_ss(2) + Delta_u_id; 

% The perturbation vector Delta_u_F2
Delta_u_F2 = u_F2_step - u_ss;

% Output Perturbation Delta_z_F2
Z_F2 = Z_det(:, :, k_step_id, 2);
Delta_z_F2 = Z_F2 - z_ss;

% --- 2. Create iddata Objects (Without InputOffset) ---

io_data_F1 = iddata(Delta_z_F1', Delta_u_F1', dt, 'TStart', t(1));
io_data_F2 = iddata(Delta_z_F2', Delta_u_F2', dt, 'TStart', t(1));

% 4b. Identify the Transfer Function Matrix G(s)
try
    % Identify G11(s) and G21(s) (outputs 1&2 from input 1)
    G_col1_tf = tfest(io_data_F1, [1 1; 1 1], [1 1; 1 1], 'Ts', 0);
    
    % Identify G12(s) and G22(s) (outputs 1&2 from input 2)
    G_col2_tf = tfest(io_data_F2, [1 1; 1 1], [1 1; 1 1], 'Ts', 0);
    
    % Assemble the final G(s) matrix
    G_tf = [G_col1_tf(:, 1), G_col2_tf(:, 2)];
    
    % Extract the gains for the display section (Section 4)
    K_final = dcgain(G_tf);
    K11_calc = K_final(1, 1);
    K12_calc = K_final(1, 2);
    K21_calc = K_final(2, 1);
    K22_calc = K_final(2, 2);

    tau11 = G_tf(1, 1).Denominator{1}(1);
    tau12 = G_tf(1, 2).Denominator{1}(1);
    tau21 = G_tf(2, 1).Denominator{1}(1);
    tau22 = G_tf(2, 2).Denominator{1}(1);
    
    disp('4. Identified Transfer Function Matrix G(s) using tfest:');
    disp(G_tf);

catch ME
    disp('tfest failed or System Identification Toolbox is not available.');
    disp('Please ensure your K and tau calculations are robust.');
    
    Z_norm_id = squeeze(Z_norm_det(:, :, k_step_id, :));
    K_ij = mean(Z_norm_id(:, round(0.9*N):N, :), 2);
    K_ij = squeeze(K_ij);
    K11_calc = K_ij(1, 1);
    K12_calc = K_ij(1, 2);
    K21_calc = K_ij(2, 1);
    K22_calc = K_ij(2, 2);

    % By inspection, assuming zero dead time (theta)
    theta11 = 0; theta12 = 0; theta21 = 0; theta22 = 0;

    % --- Direct Paths (G11, G22) - 63.2% Method ---
    
    % Tau 11
    Target_11 = K11_calc * 0.632;
    k_63_11 = find(Z_norm_id(1, :, 1) >= Target_11, 1, 'first'); 
    
    if isempty(k_63_11)
        warning('Target 11 (63.2%% of K11) not reached. Using end time for tau11. Ensure N is large.');
        t_63_11 = t(end);
    else
        t_63_11 = t(k_63_11);
    end
    tau11 = t_63_11 - t_step - theta11;
    
    % Tau 22
    Target_22 = K22_calc * 0.632;
    k_63_22 = find(Z_norm_id(2, :, 2) >= Target_22, 1, 'first');
    
    if isempty(k_63_22)
        warning('Target 22 (63.2%% of K22) not reached. Using end time for tau22. Ensure N is large.');
        t_63_22 = t(end);
    else
        t_63_22 = t(k_63_22);
    end
    tau22 = t_63_22 - t_step - theta22;
    
    % --- Cross-Coupling Paths (G12, G21) - Tangential Slope Method ---
    % Method: tau = K / S, where S is the maximum slope (rate of change)
    
    % Tau 12
    Response_12 = Z_norm_id(1, k_step:N, 2); 
    Derivative_12 = diff(Response_12) / dt;
    S_12 = max(abs(Derivative_12));
    if S_12 == 0
        warning('Max slope S_12 is zero. Tau12 is undefined. Setting to T_end.');
        tau12 = t(end); 
    else
        tau12 = abs(K12_calc) / S_12; 
    end

    % Tau 21
    Response_21 = Z_norm_id(2, k_step:N, 1);
    Derivative_21 = diff(Response_21) / dt;
    S_21 = max(abs(Derivative_21));
    if S_21 == 0
        warning('Max slope S_21 is zero. Tau21 is undefined. Setting to T_end.');
        tau21 = t(end); 
    else
        tau21 = abs(K21_calc) / S_21;
    end

    disp('Tau parameters:')
    disp([tau11 tau12 tau21 tau22])
end


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

k_start = k_step; 
k_end = N;
N_M = k_end - k_start + 1; % Length of the Markov Parameter sequence

% Use the deterministic, 10% step normalized response (Z_norm_id) for M_k.
% M now stores only the sequence AFTER the step (M[output_i, time, input_j])
M = zeros(ny, N_M, nu); 
M_t = t(k_start:k_end);

for i = 1:nu
    for l = 1:ny
        % Get the normalized response *only after the step*
        Response_segment = Z_norm_id(l, k_start:k_end, i);
        M(l, 1:N_M-1, i) = diff(Response_segment) / dt;
    end
end

figure();
sgtitle('Impulse Response Coefficients (Markov Parameters)');

% The time vector must be shifted to start at 0
t_shift = (M_t(2:N_M) - M_t(1))/60; % Time in minutes, relative to the step

for i = 1:nu
    for l = 1:ny
        subplot(nu, ny, (i-1)*ny + l);
        hold on;
        
        % Plot starts from t=dt (i.e., the first change after the impulse)
        plot(t_shift, M(l, 1:N_M-1, i), 'DisplayName', 'Markov Parameters');
        
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