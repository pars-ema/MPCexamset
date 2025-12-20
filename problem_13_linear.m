%   Problem based on Problem 10 implementation (input-constrainted MPC + KF)
%   Aim is to include **economic constraints** on defined problem

%%  1. Problem formulation

%{
1. OBJECTIVE FUNCTIONS (Bi-criterion Structure)
   The problem follows a Bi-criterion optimization structure:
   phi = alpha * phi_eco + (1 - alpha) * phi_reg, where alpha in [0, 1]
   
   Where:
   - phi_reg: Least Squares Objective (Tracking)
     phi_reg = 0.5 * sum(||x_k - x_bar||_Q^2 + ||u_k - u_bar||_R^2)
   - phi_eco: Economic Objective (cost, profit, ...)
     phi_eco = sum(l_k(x_k, u_k, theta))

However, the approach chosen is another one based on the gradient of the
given function. The given alternative is explained below.

2. GRADIENT-BASED IMPLEMENTATION
   The economic cost is integrated via the linear gradient of phi_eco. 
   Given the pumping cost rate P(t) = c1*F1 + c2*F2, the discrete 
   economic gradient is derived as:
   
   f_eco = grad(phi_eco) = kron(ones(N,1), [c1; c2] * Ts)

3. TOTAL GRADIENT ASSEMBLY
   Total Gradient (g_U) = (1 - beta) * grad(phi_reg) + beta * f_eco
   
   *Note: beta in the code corresponds to alpha in the bi-criterion 
   formulation.
%}


clear;clc;close all;

%   Run problem 10 (linear MPC) to obtain all specifications.
%   Then, close plots with 'close all'
run("problem_10_offset_free.m")
close all;


%% ========================================================================
%   Economic MPC - linear KF discrete cost function & eco. parameters
% ========================================================================

%   Cost function (Economic MPC) - Definition of cost for (F1,F2) flow
%   rates, also based on the flow time window for discretization.
cost_linear = [1;0.5];
%cost_ss_term = (cost_linear(1)*u_s(1) + cost_linear(2)*u_s(2)) * Ts;
%cost_linear_func = @(u,dt) cost_ss_term + (cost_linear(1)*u(1) + ...
%                    cost_linear(2)*u(2)) * dt;

%   This coefficient alpha determines the reg-eco terms' relationship in 
%   the J_eco economical cost function calculation.
%alpha = 0.5;

%   Coefficient for economic cost addition in J. The higher the beta coefficient,
%   the more influence the economic cost term will have in the cost minimization.
beta = 0.3;


%% ========================================================================
%   Closed-loop simulation (OFFSET-FREE MPC + static augmented KF)
% ========================================================================

sim = simulateMPC_Problem13( ...
    A,B,E,C, ...
    A_e,B_e,C_e, ...
    K_e, ...
    Phi_z_x,Phi_z_u, ...
    Wz_big,Wu_big,Wdu_big, ...
    E_du, ...
    u_s, ...
    u_min_stack,u_max_stack, ...
    du_min_stack,du_max_stack, ...
    Z_max_stack, ...
    z_ref_abs,z_ref, ...
    z_min_abs,z_max_abs, ...
    zs, ...
    d_dev, ...
    Tsim,Ts,R_meas,cost_linear,beta);


%% ========================================================================
%   Plots – outputs, inputs, constraints, disturbances
% ========================================================================

plotResultsProblem13(t,sim, ...
    z_ref_abs,z_ref, ...
    u_min_abs,u_max_abs, ...
    du_min_abs,du_max_abs, ...
    z_min_abs,z_max_abs, ...
    k_step, d_step_val);


%% ========================================================================
% ======================== LOCAL FUNCTIONS ================================
% ========================================================================

function sim = simulateMPC_Problem13( ...
    A,B,E,C, ...
    A_e,B_e,C_e, ...
    K_e, ...
    Phi_z_x,Phi_z_u, ...
    Wz_big,Wu_big,Wdu_big, ...
    E_du, ...
    u_s, ...
    u_min_stack,u_max_stack, ...
    du_min_stack,du_max_stack, ...
    Z_max_stack, ...
    z_ref_abs,z_ref, ...
    z_min_abs,z_max_abs, ...
    zs, ...
    d_dev, ...
    Tsim,Ts,R_meas,cost_linear,beta)

    [nx,nu] = size(B);
    nd = size(E,2);
    ny = size(C,1);
    nx_e = nx + nd;

    N  = size(Phi_z_u,1) / ny;   % prediction horizon
    nU = N*nu;
    nZ = N*ny;
    nEta = nZ;                   % one slack per predicted output

    % Slack penalty (softness) – big = stricter
    w_eta = 1e3;
    Weta_big = w_eta * eye(nEta);

    % Precompute constant Hessian part for U
    H_z  = Phi_z_u' * Wz_big  * Phi_z_u;
    H_u  = Wu_big;
    H_du = E_du'   * Wdu_big * E_du;
    H_UU = H_z + H_u + H_du;

    % Box bounds for η (slack)
    eta_max = 50;                    % [cm] extra allowed above z_max_dev
    eta_min_vec = zeros(nEta,1);
    eta_max_vec = eta_max * ones(nEta,1);

    % Plant noise (linear)
    Qplant = 1e-4 * eye(nx);
    sigma_meas = sqrt(diag(R_meas));

    % Deviation bounds for outputs
    z_min_dev = z_min_abs - zs;
    z_max_dev = z_max_abs - zs;
    Z_max_stack_dev = Z_max_stack; %#ok<NASGU> (name clarity)

    % Initialise
    x_true = zeros(nx,1);            % deviation state
    d_true = zeros(nd,1);            % deviation disturbance
    x_hat_e = zeros(nx_e,1);         % augmented estimate [x_hat; d_hat]
    u_prev = zeros(nu,1);

    % Logs
    z_dev = zeros(ny,Tsim);
    z_abs = zeros(ny,Tsim);
    u_dev = zeros(nu,Tsim);
    u_abs = zeros(nu,Tsim);
    d_hat = zeros(nd,Tsim);
    eps0  = zeros(ny,Tsim);          % first-step slack

    nu_total = nU + nEta;

    for k = 1:Tsim

        % Disturbance for this step (deviation)
        d_true = d_dev(:,k);

        % ----- True plant update -----
        if k > 1
            w = mvnrnd(zeros(nx,1),Qplant)';
            x_true = A*x_true + B*u_prev + E*d_dev(:,k-1) + w;
        end

        z_dev(:,k) = C*x_true;
        z_abs(:,k) = zs + z_dev(:,k);

        % Measurement with noise (deviation)
        v = sigma_meas .* randn(ny,1);
        y_meas = z_dev(:,k) + v;

        % ----- Static augmented Kalman filter (offset-free) -----
        % Augmented measurement model: y = C_e x_e, with C_e = [C  0]
        y_pred = C_e*x_hat_e;
        innov  = y_meas - y_pred;
        x_hat_e = x_hat_e + K_e*innov;

        % Store estimated disturbances
        d_hat(:,k) = x_hat_e(nx+1:end);

        % ----- Build horizon references (deviation) -----
        if k+N-1 <= Tsim
            zbar = z_ref(:,k:k+N-1);
            ubar = zeros(nu,N);          % nominal input deviations = 0
        else
            last = Tsim-k+1;
            zbar = [z_ref(:,k:end), repmat(z_ref(:,end),1,N-last)];
            ubar = zeros(nu,N);
        end

        Zbar = zbar(:);
        Ubar = ubar(:);

        % ----- Cost gradient for U -----
        y0 = Phi_z_x*x_hat_e;          % predicted output dev offset term

        f_z  = Phi_z_u' * Wz_big  * (y0 - Zbar);
        f_u  = -Wu_big * Ubar;

        Ed_prev = [u_prev; zeros((N-1)*nu,1)];
        f_du = -E_du' * Wdu_big * Ed_prev;

        %---------------------------------------------------------
        %   Economic approach inclusion (f_eco economic gradient)
        f_eco = kron(ones(N,1),cost_linear(:)*Ts);

        %   Multi-objective Gradient Assembly
        %   beta = 1 means purely economic, beta = 0 means purely tracking
        g_U = (1-beta)*(f_z + f_u + f_du) + beta*(f_eco);
        %---------------------------------------------------------

        % Full decision vector w = [U; eta]
        H = blkdiag(H_UU, Weta_big);
        f = [g_U; zeros(nEta,1)];

        % ----- Box constraints: U and eta -----
        l = [u_min_stack; eta_min_vec];
        u = [u_max_stack; eta_max_vec];

        % ----- Move constraints: ΔUmin <= E_du U - Ed_prev <= ΔUmax -----
        A1  = [E_du, zeros(N*nu,nEta)];
        bl1 = du_min_stack + Ed_prev;
        bu1 = du_max_stack + Ed_prev;

        % ----- Soft output upper bounds: z_pred_dev <= z_max_dev + eta -----
        % z_pred_dev = y0 + Phi_z_u U
        %   y0 + Phi_z_u U <= Z_max_stack + eta
        % => Phi_z_u U - I*eta <= Z_max_stack - y0
        Zmax_stack = Z_max_stack;          % (N*ny x 1) dev-bounds-stack
        b_soft = Zmax_stack - y0;

        bigM = 1e6;
        A2  = [Phi_z_u, -eye(nZ)];
        bl2 = -bigM * ones(nZ,1);         % effectively one-sided
        bu2 = b_soft;

        % Stack inequalities
        Aqp   = [A1; A2];
        bl_ineq = [bl1; bl2];
        bu_ineq = [bu1; bu2];

        % ----- Solve QP with problem_7 -----
        xinit = zeros(nu_total,1);
        [w_opt,~] = problem_7(H,f,l,u,Aqp,bl_ineq,bu_ineq,xinit);

        Uopt   = w_opt(1:nU);
        etaopt = w_opt(nU+1:end);         %#ok<NASGU> (full seq, only first shown)

        u_k = Uopt(1:nu);

        u_dev(:,k) = u_k;
        u_abs(:,k) = u_s + u_k;

        % First-step slack for each output (for plotting)
        eps0(:,k) = etaopt(1:ny);

        % Time update of KF
        x_hat_e = A_e*x_hat_e + B_e*u_k;

        u_prev = u_k;
    end

    % Pack results
    sim.z_dev = z_dev;
    sim.z_abs = z_abs;
    sim.u_dev = u_dev;
    sim.u_abs = u_abs;
    sim.d_hat = d_hat;
    sim.d_true = d_dev;
    sim.eps = eps0;
end

% -------------------------------------------------------------------------
function plotResultsProblem13(t,sim, ...
    z_ref_abs,z_ref, ...
    u_min_abs,u_max_abs, ...
    du_min_abs,du_max_abs, ...
    z_min_abs,z_max_abs, ...
    k_step, d_step_val)

    z_abs = sim.z_abs;
    z_dev = sim.z_dev;
    u_abs = sim.u_abs;
    u_dev = sim.u_dev;
    d_hat = sim.d_hat;
    d_true = sim.d_true;
    eps0  = sim.eps;

    %% Absolute outputs with soft bounds
    figure('Name','P10 – Absolute Outputs (Offset-free, Soft Constraints)');
    subplot(2,1,1); hold on; grid on;
    plot(t, z_abs(1,:), 'b','LineWidth',1.5);
    plot(t, z_ref_abs(1,:), 'k--','LineWidth',1.2);
    yline(z_min_abs(1),'g:','LineWidth',1.2);
    yline(z_max_abs(1),'g:','LineWidth',1.2);
    xline(t(k_step),'r--','d-step','LabelVerticalAlignment','bottom');
    ylabel('h_1 [cm]');
    legend('h_1','h_{1,ref}','h_{1,min}^{soft}','h_{1,max}^{soft}','Location','best');
    title('Tank 1 – Absolute Height with Soft Bounds');

    subplot(2,1,2); hold on; grid on;
    plot(t, z_abs(2,:), 'b','LineWidth',1.5);
    plot(t, z_ref_abs(2,:), 'k--','LineWidth',1.2);
    yline(z_min_abs(2),'g:','LineWidth',1.2);
    yline(z_max_abs(2),'g:','LineWidth',1.2);
    xline(t(k_step),'r--','d-step','LabelVerticalAlignment','bottom');
    ylabel('h_2 [cm]'); xlabel('Time [s]');
    legend('h_2','h_{2,ref}','h_{2,min}^{soft}','h_{2,max}^{soft}','Location','best');
    title('Tank 2 – Absolute Height with Soft Bounds');

    %% Deviation outputs
    figure('Name','P10 – Deviation Outputs (Offset-free)');
    subplot(2,1,1); hold on; grid on;
    plot(t, z_dev(1,:), 'b','LineWidth',1.5);
    plot(t, z_ref(1,:), 'k--','LineWidth',1.2);
    xline(t(k_step),'r--','d-step');
    ylabel('z_1 dev [cm]');
    legend('z_1','z_{1,ref}','Location','best');
    title('Output 1 (deviation)');

    subplot(2,1,2); hold on; grid on;
    plot(t, z_dev(2,:), 'b','LineWidth',1.5);
    plot(t, z_ref(2,:), 'k--','LineWidth',1.2);
    xline(t(k_step),'r--','d-step');
    ylabel('z_2 dev [cm]'); xlabel('Time [s]');
    legend('z_2','z_{2,ref}','Location','best');
    title('Output 2 (deviation)');

    %% Absolute inputs + hard bounds
    figure('Name','P10 – Absolute Inputs (Offset-free)');
    hold on; grid on;
    plot(t, u_abs(1,:), 'r','LineWidth',1.5);
    plot(t, u_abs(2,:), 'm--','LineWidth',1.5);

    yline(u_min_abs(1),'b--','LineWidth',1.2);
    yline(u_max_abs(1),'b--','LineWidth',1.2);
    yline(u_min_abs(2),'g--','LineWidth',1.2);
    yline(u_max_abs(2),'g--','LineWidth',1.2);
    xline(t(k_step),'r--','d-step');

    xlabel('Time [s]');
    ylabel('Flow [cm^3/s]');
    legend('u_1','u_2', ...
           'u_{1,min}','u_{1,max}', ...
           'u_{2,min}','u_{2,max}','Location','best');
    title('Absolute Pump Flows with Input Constraints');

    %% Deviation inputs
    figure('Name','P10 – Deviation Inputs (Offset-free)');
    hold on; grid on;
    stairs(t, u_dev(1,:), 'r','LineWidth',1.5);
    stairs(t, u_dev(2,:), 'm--','LineWidth',1.5);
    xline(t(k_step),'r--','d-step');
    xlabel('Time [s]');
    ylabel('\Delta u [cm^3/s]');
    legend('\Delta u_1','\Delta u_2','Location','best');
    title('Deviation Inputs');

    %% Input moves + move constraints
    du1 = [u_dev(1,1), diff(u_dev(1,:))];
    du2 = [u_dev(2,1), diff(u_dev(2,:))];

    figure('Name','P10 – Input Moves (Offset-free)');
    hold on; grid on;
    stairs(t, du1, 'r','LineWidth',1.5);
    stairs(t, du2, 'm--','LineWidth',1.5);

    yline(du_max_abs(1),'b--');
    yline(du_min_abs(1),'b--');
    yline(du_max_abs(2),'g--');
    yline(du_min_abs(2),'g--');
    xline(t(k_step),'r--','d-step');

    xlabel('Time [s]');
    ylabel('\Delta u [cm^3/s]');
    legend('\Delta u_1','\Delta u_2', ...
           '\Delta u_{1,max}','\Delta u_{1,min}', ...
           '\Delta u_{2,max}','\Delta u_{2,min}','Location','best');
    title('Input Moves and Move Constraints');

    %% Disturbances true vs estimated (offset-free behaviour)
    figure('Name','P10 – Disturbances (Offset-free)');
    subplot(2,1,1); hold on; grid on;
    plot(t, d_true(1,:), 'k--','LineWidth',1.4,'DisplayName','True \DeltaF_3');
    plot(t, d_hat(1,:),  'b',  'LineWidth',1.5,'DisplayName','\hat{\DeltaF}_3');
    xline(t(k_step),'r--','d-step');
    ylabel('\DeltaF_3 [cm^3/s]');
    legend('Location','best');
    title(sprintf('F_3 disturbance step = %.1f cm^3/s', d_step_val(1)));

    subplot(2,1,2); hold on; grid on;
    plot(t, d_true(2,:), 'k--','LineWidth',1.4,'DisplayName','True \DeltaF_4');
    plot(t, d_hat(2,:),  'b',  'LineWidth',1.5,'DisplayName','\hat{\DeltaF}_4');
    xline(t(k_step),'r--','d-step');
    ylabel('\DeltaF_4 [cm^3/s]'); xlabel('Time [s]');
    legend('Location','best');
    title('F_4 disturbance (here no step)');
end