clc; clear; close all
% -----------------------------------------------------------
% Parameters
% -----------------------------------------------------------
a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;        % [cm2]
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327; % [cm2]
g = 981; rho = 1.00;
gamma1 = 0.58; gamma2 = 0.72;
p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2; rho];

F1 = 300; F2 = 300;
% Mean disturbances
F3 = 250; F4 = 250;           % [cm3/s]
d = [F3; F4];

% Steady-state inputs (nominal pump flows)
u_s = [300; 300];             % [cm3/s]

% ------------------------------------------------------------
% Compute steady state automatically using fsolve
% ------------------------------------------------------------
xs0 = 500*ones(4,1);  % initial guess
opts = optimoptions('fsolve','Display','none','FunctionTolerance',1e-12,'StepTolerance',1e-12);
xs = fsolve(@(x) FourTankSystemSteadyResidual(x,u_s,d,p), xs0, opts);

% Compute steady-state heights for reference
hs = xs ./ (rho * [A1; A2; A3; A4]);

disp('Steady-state heights [cm]:'); disp(hs.');
disp('Steady-state masses [g]:'); disp(xs.');


%% CONTINUOUS-TIME TRANSFER FUNCTIONS FROM THE LINEARIZED MODEL (mass-based)

% 0) Unpack parameters (matching your vector p)
a1=p(1); a2=p(2); a3=p(3); a4=p(4);
A1=p(5); A2=p(6); A3=p(7); A4=p(8);
g =p(9); gamma1=p(10); gamma2=p(11); rho=p(12);

% 1) Compute linearization coefficients at the steady state xs (masses)
% Torricelli constants
c1 = a1*sqrt(2*g/(rho*A1));
c2 = a2*sqrt(2*g/(rho*A2));
c3 = a3*sqrt(2*g/(rho*A3));
c4 = a4*sqrt(2*g/(rho*A4));

beta1 = c1/(2*sqrt(xs(1)));
beta2 = c2/(2*sqrt(xs(2)));
beta3 = c3/(2*sqrt(xs(3)));
beta4 = c4/(2*sqrt(xs(4)));


%% Build Jacobian-based state-space matrices (mass-based linearization) described in the report

% The Jacobian of f(x,u,d,p) with respect to states (evaluated at steady state)
A = rho * [
   -beta1, 0,       beta3,   0;
        0, -beta2,       0,  beta4;
        0,      0,  -beta3,   0;
        0,      0,       0,  -beta4
];

% Jacobian with respect to inputs (pumps F1, F2)
B = rho * [
    gamma1,        0;
         0,   gamma2;
         0, 1-gamma2;
   1-gamma1,       0
];

% Jacobian with respect to disturbances (F3, F4)
E = rho * [
    0, 0;
    0, 0;
    1, 0;
    0, 1
];

% Output matrix: maps masses to measured levels y = [h1; h2]
C = [
    1/(rho*A1), 0, 0, 0;
    0, 1/(rho*A2), 0, 0
];



%%chatGPT proposed following so i dont really understand it. I would use
%%classical:
% figure
% sys_u = ss(A,B,C,zeros(2,2)); % From inputs F1,F2 to outputs h1,h2
% Gyu = tf(sys_u); % Transfer functions G_ij(s)
% sys_d = ss(A,E,C,zeros(2,2)); % From disturbances F3,F4 to outputs
% Gyd = tf(sys_d);
% step(Gyu)

a1_lin = rho*beta1;  a2_lin = rho*beta2;  a3_lin = rho*beta3;  a4_lin = rho*beta4;

% Time constants
tau1 = 1/a1_lin;  tau2 = 1/a2_lin;  tau3 = 1/a3_lin;  tau4 = 1/a4_lin;

% DC gains in canonical form
K11 = gamma1        /(A1*a1_lin);
K12 = (1-gamma2)    /(A1*a1_lin);
K21 = (1-gamma1)    /(A2*a2_lin);
K22 = gamma2        /(A2*a2_lin);

Kd1 = 1/(A1*a1_lin);
Kd2 = 1/(A2*a2_lin);

% Build G(s) = K/(tau s + 1) or K/((tau1 s+1)(tau2 s+1))
s = tf('s');
G11 = K11 / (tau1*s + 1);
G12 = K12 / ((tau1*s + 1)*(tau3*s + 1));
G21 = K21 / ((tau2*s + 1)*(tau4*s + 1));
G22 = K22 / (tau2*s + 1);

Gyu_canonical = [G11 G12; G21 G22];

% Disturbance channels (optional)
Gy1d1 = Kd1 / ((tau1*s + 1)*(tau3*s + 1));
Gy2d2 = Kd2 / ((tau2*s + 1)*(tau4*s + 1));
Gyd_canonical = [Gy1d1 0; 0 Gy2d2];

% Pretty print
disp('Canonical Gyu(s):'); Gyu_canonical
disp('Canonical Gyd(s):'); Gyd_canonical
% 5) (Optional) plot step responses for each SISO element of Gyu
t0 = 0.0;        % [s]
t_final = 60*30;      % [s]
Ts = 1;          % [s]
t = t0:Ts:t_final;
figure('Name','Step responses of G_{yu}(s)'); 
subplot(2,2,1); step(Gyu_canonical(1,1),t); grid on; title('G_{11}: h_1 \leftarrow F_1','Interpreter','tex', LineWidth=1.5);
ylim([0 0.5]);
subplot(2,2,2); step(Gyu_canonical(1,2),t); grid on; title('G_{12}: h_1 \leftarrow F_2','Interpreter','tex', LineWidth=1.5);
ylim([0 0.5]);
subplot(2,2,3); step(Gyu_canonical(2,1),t); grid on; title('G_{21}: h_2 \leftarrow F_1','Interpreter','tex', LineWidth=1.5);
ylim([0 0.5]);
subplot(2,2,4); step(Gyu_canonical(2,2),t); grid on; title('G_{22}: h_2 \leftarrow F_2','Interpreter','tex', LineWidth=1.5);
ylim([0 0.5]);
outputFolder = fullfile('figures', 'problem_5', 'step_response_linearized','deterministic'); 
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
exportgraphics(gcf, fullfile(outputFolder, 'step_response.pdf'), 'ContentType', 'vector');

%%

figure
sys_u = ss(A,B,C,zeros(2,2)); % From inputs F1,F2 to outputs h1,h2
Gyu = tf(sys_u); % Transfer functions G_ij(s)
% sys_d = ss(A,E,C,zeros(2,2)); % From disturbances F3,F4 to outputs
% Gyd = tf(sys_d);
step(Gyu)








%% dicretization with mimoctf2d
%TO DO
% ------------------------------------------------------------
%  Discretization of the linearized 4-tank model with mimoctf2dss
% ------------------------------------------------------------
Ts   = 10;        % sampling time [s]
Nmax = 100;      % maximum state dimension
tol  = 1e-8;     
lambda = zeros(2,2); % assume no pure time delay

% Extract numerators / denominators from your continuous model Gyu
num = cell(2,2);
den = cell(2,2);
for i = 1:2
    for j = 1:2
        [num{i,j}, den{i,j}] = tfdata(Gyu(i,j),'v');
    end
end

% Convert continuous-time MIMO TF to discrete-time state-space
[Ad,Bd,Cd,Dd,sH] = mimoctf2dss(num,den,lambda,Ts,Nmax,tol);

% Display results
disp('Discrete-time state-space matrices:');
Ad, Bd, Cd, Dd

% Save them to a .mat file
save('discrete_model_from_linearization.mat', 'Ad', 'Bd', 'Cd', 'Dd', 'sH');


% Optional: plot Hankel singular values (to inspect order content)
figure; semilogy(sH,'b.-','LineWidth',2,'MarkerSize',20);
grid on;
xlabel('State number'); ylabel('Singular value');
title('Hankel Singular Values of Discrete 4-Tank Model');

Nimp = 1800/Ts;
ny = size(Cd,1);
nu = size(Bd,2);

H_det = zeros(ny, nu, Nimp);
H_det(:,:,1) = Dd;
for k = 2:Nimp
    H_det(:,:,k) = Cd * (Ad^(k-2)) * Bd;
end

% Save deterministic (exact) Markov parameters
save('Markov_parameters_problem_5.mat','H_det','Ts');
disp('✅ Saved deterministic Markov parameters to Markov_parameters_problem_5.mat');

%%comparison with parameters from the previous problem
%% Compare Markov parameters: Identified vs Deterministic
load('H_identified_problem_4.mat','H_identified_problem_4','Ts');
load('Markov_parameters_problem_5.mat','H_det');
%%
Nimp = min(size(H_identified_problem_4,3), size(H_det,3));
k = 0:Nimp-1;

figure('Name','Comparison of Markov Parameters','Position',[100 100 1200 800]);
titles = {'$H_{11}$','$H_{12}$','$H_{21}$','$H_{22}$'};
for i = 1:2
    for j = 1:2
        subplot(2,2,(i-1)*2+j); hold on; grid on;
        plot(k, squeeze(H_det(i,j,1:Nimp)), 'bo', 'LineWidth',1.5, 'DisplayName','Deterministic linearized');
        plot(k, squeeze(H_identified_problem_4(i,j,1:Nimp)), 'rx', 'LineWidth',1.3, 'DisplayName','Identified TF');
        xlabel('Sample k');
        ylabel('Amplitude');
        title(titles{(i-1)*2+j}, 'Interpreter','latex', 'FontSize',14);
        legend('Location','best');
    end
end
sgtitle('Comparison of Markov Parameters: Deterministic vs Identified Model', ...
        'Interpreter','latex','FontSize',16);

outputFolder = fullfile('figures','problem_6','Markov_comparison');
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
exportgraphics(gcf, fullfile(outputFolder,'Markov_comparison.pdf'),'ContentType','vector');
disp('✅ Markov parameter comparison figure saved.');



%% Tranfer function in nice form
% Convert from SS to TF

% Extract each transfer function
G11 = Gyu(1,1);
G12 = Gyu(1,2);
G21 = Gyu(2,1);
G22 = Gyu(2,2);



% Compute for all G_ij
[K11, tau11] = extractParams(G11);
[K12, tau12] = extractParams(G12);
[K21, tau21] = extractParams(G21);
[K22, tau22] = extractParams(G22);

% Print in canonical form
fprintf('\n--- Canonical Transfer Function Forms from linearization ---\n');
fprintf('G11(s) = %.4e / (%.2fs + 1)\n', K11, tau11);
fprintf('G12(s) = %.4e / ((%.2fs + 1)(%.2fs + 1))\n', K12, tau12(1), tau12(2));
fprintf('G21(s) = %.4e / ((%.2fs + 1)(%.2fs + 1))\n', K21, tau21(1), tau21(2));
fprintf('G22(s) = %.4e / (%.2fs + 1)\n', K22, tau22);

%%

% ------------------------------------------------------------
% Supporting residual function
% ------------------------------------------------------------
function r = FourTankSystemSteadyResidual(x, u, d, p)
r = FourTankSystemModified(0, x, u, d, p); % residual f(x)=0
end

% Helper to extract K, tau(s)
function [K, tau] = extractParams(G)
    [num, den] = tfdata(G, 'v');
    % Normalize denominator
    den = den / den(end); 
    poles = roots(den);
    taus = -1 ./ real(poles);
    K = dcgain(G);
    tau = sort(taus);
end



