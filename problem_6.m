clear; clc; close all;

% --- Load discrete-time models from your previous work ---
load('discrete_model_from_step.mat',          'Ad','Bd','Cd','Dd');  Ad_step=Ad; Bd_step=Bd; Cd_step=Cd; Dd_step=Dd;
load('discrete_model_from_linearization.mat', 'Ad','Bd','Cd','Dd');  Ad_lin =Ad; Bd_lin =Bd; Cd_lin =Cd; Dd_lin =Dd;

Ts = 10;                        % <- set your sampling time if needed
nxS = size(Ad_step,1); nyS = size(Cd_step,1); nuS = size(Bd_step,2);
nxL = size(Ad_lin ,1); nyL = size(Cd_lin ,1); nuL = size(Bd_lin ,2);
N   = 800;                     % simulation length


%%
% --- Tuning (start conservative, scale later) ---
Q_step = 12.5*eye(nxS);   R_step = 2*eye(nyS);  S_step = zeros(nxS,nyS);
Q_lin  = 1e-6*eye(nxL);   R_lin  = 1e-4*eye(nyL);  S_lin  = zeros(nxL,nyL);

% --- Static (steady-state) Kalman filter gains (dlqe) ---
% Identified model
[Lss_step, Pss_step, ~] = dlqe(Ad_step, eye(nxS), Cd_step, Q_step, R_step, S_step);
% Linearized model
[Lss_lin , Pss_lin , ~] = dlqe(Ad_lin , eye(nxL), Cd_lin , Q_lin ,  R_lin ,  S_lin );

% --- Helper: dynamic KF recursion ---
kf_dyn = @(A,B,C,Q,R,u,y,x0,P0) ...
    local_dynamic_kf(A,B,C,Q,R,u,y,x0,P0);

% --- Helper: static KF recursion (fixed L) ---
kf_ss  = @(A,B,C,L,u,y,x0) ...
    local_static_kf(A,B,C,L,u,y,x0);


% --- Choose which model generates the "plant" for this experiment ---
Atrue = Ad_step; Btrue = Bd_step; Ctrue = Cd_step; Dtrue = Dd_step;
nx = size(Atrue,1); ny = size(Ctrue,1); nu = size(Btrue,2);

% Inputs (you can replace with your real inputs)
u = zeros(nu,N);
u(1,1:N) = 250; u(2,1:N) = 250;   % small excitation

% Generate a baseline stochastic scenario (no step changes)
Qtrue = 1e-7*eye(nx);   Rtrue = 1e-4*eye(ny);
xtrue = zeros(nx,N); ymeas = zeros(ny,N);
rng(1);
for k=2:N
    xtrue(:,k) = Atrue*xtrue(:,k-1) + Btrue*u(:,k-1) + chol(Qtrue,'lower')*randn(nx,1);
    ymeas(:,k) = Ctrue*xtrue(:,k)   + Dtrue*u(:,k)   + chol(Rtrue,'lower')*randn(ny,1);
end

% Run dynamic KF (time-varying)
x0 = zeros(nx,1); P0 = 1e-2*eye(nx);
[xhat_dyn, ~] = kf_dyn(Atrue,Btrue,Ctrue,Q_step,R_step,u,ymeas,x0,P0);

% Run static KF (steady-state)
xhat_ss = kf_ss(Atrue,Btrue,Ctrue,Lss_step,u,ymeas,x0);

% Compare
figure; t=0:N-1;
subplot(2,1,1); plot(t,xtrue(1,:),'-', t,xhat_dyn(1,:),'--', t,xhat_ss(1,:),':','LineWidth',1.25);
xlabel('k'); ylabel('x_1'); legend('True','Dynamic KF','Static KF'); grid on;
subplot(2,1,2); plot(t,xtrue(2,:),'-', t,xhat_dyn(2,:),'--', t,xhat_ss(2,:),':','LineWidth',1.25);
xlabel('k'); ylabel('x_2'); legend('True','Dynamic KF','Static KF'); grid on;
title('Baseline noise (no step disturbances)');


%%
% --- Augmentation for output-bias random walk ---
% x_aug = [x; d],  y = [I]*[C  I]*[x; d] + v
Abar = [Atrue, zeros(nx,ny); zeros(ny,nx), eye(ny)];
Bbar = [Btrue; zeros(ny,nu)];
Cbar = [Ctrue, eye(ny)];
Dbar = [Dtrue];

% Noise covariances: process noise includes both x and bias d
Qx   = Q_step;              % state process noise (tune)
Qd   = (1e-7)*eye(ny);      % bias random walk (small to allow slow steps)
Qbar = blkdiag(Qx, Qd);
Rbar = R_step;

% Initials
x0bar = zeros(nx+ny,1);  P0bar = blkdiag(1e-2*eye(nx), 1e-2*eye(ny));

% Build a test with actual step change in output bias at k=350
y_meas_step = ymeas;
y_meas_step(:,350:end) = y_meas_step(:,350:end) + [0.3; -0.25];  % step biases

% Dynamic KF on augmented system
[xhat_dyn_bar, ~] = kf_dyn(Abar,Bbar,Cbar,Qbar,Rbar,u,y_meas_step,x0bar,P0bar);

% Static (steady-state) KF on augmented system
[Lss_bar, ~, ~] = dlqe(Abar, eye(nx+ny), Cbar, Qbar, Rbar, zeros(nx+ny,ny));
xhat_ss_bar = kf_ss(Abar,Bbar,Cbar,Lss_bar,u,y_meas_step,x0bar);

% Extract state estimates (first nx are the physical states; last ny are biases)
xhat_dyn_of = xhat_dyn_bar(1:nx, :);
dbias_dyn   = xhat_dyn_bar(nx+1:end, :);
xhat_ss_of  = xhat_ss_bar(1:nx, :);
dbias_ss    = xhat_ss_bar(nx+1:end, :);

% Compare
figure; t=0:N-1;
subplot(2,1,1); plot(t,xtrue(1,:),'-', t,xhat_dyn_of(1,:),'--', t,xhat_ss_of(1,:),':','LineWidth',1.25);
xlabel('k'); ylabel('x_1'); legend('True','Dyn OF-KF','SS OF-KF'); grid on;
subplot(2,1,2); plot(t, y_meas_step(1,:) - (Ctrue(1,:)*xhat_dyn_of), 'LineWidth',1.25); hold on;
plot(t, dbias_dyn(1,:), '--', 'LineWidth',1.25);
legend('Residual (y-Cx)','Estimated bias d_1'); xlabel('k'); grid on;
title('Offset-free augmentation (output-bias random walk)');

%%
function [xhat, P_hist] = local_dynamic_kf(A,B,C,Q,R,u,y,x0,P0)
nx = size(A,1); N = size(y,2);
xhat = zeros(nx,N); xhat(:,1) = x0; P = P0; P_hist = zeros(nx,nx,N); P_hist(:,:,1)=P;
for k=2:N
    % Prediction
    xpred = A*xhat(:,k-1) + B*u(:,k-1);
    Ppred = A*P*A' + Q;
    % Update
    Re = C*Ppred*C' + R;
    K  = (Ppred*C')/Re;
    e  = y(:,k) - C*xpred;
    xhat(:,k) = xpred + K*e;
    P  = (eye(nx)-K*C)*Ppred;
    P_hist(:,:,k)=P;
end
end

function xhat = local_static_kf(A,B,C,L,u,y,x0)
nx = size(A,1); N = size(y,2);
xhat = zeros(nx,N); xhat(:,1)=x0;
for k=2:N
    xpred = A*xhat(:,k-1) + B*u(:,k-1);
    e     = y(:,k) - C*xpred;
    xhat(:,k) = xpred + L*e;      % fixed gain
end
end


