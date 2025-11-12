close all
clc

% to run this, you have to run either problem 4 or 5 before hand

Ts = 1;

nx = 4; nu = 2;
fhandle = @(x,u) FourTankSystemModified(0, x, u, d, p);

Acont = zeros(nx,nx);
Bcont = zeros(nx,nu);
h = 1e-6;

for i = 1:nx
    dx = zeros(nx,1); dx(i) = h;
    fp = fhandle(xs + dx, u_s);
    fm = fhandle(xs - dx, u_s);
    Acont(:,i) = (fp - fm)/(2*h);
end

for j = 1:nu
    du = zeros(nu,1); du(j) = h;
    fp = fhandle(xs, u_s + du);
    fm = fhandle(xs, u_s - du);
    Bcont(:,j) = (fp - fm)/(2*h);
end

% output: heights of all 4 tanks
S = rho * [A1; A2; A3; A4];
Ccont = diag(1./S);
Dcont = zeros(4,2);

% discretize
sysd = c2d(ss(Acont,Bcont,Ccont,Dcont), Ts);

Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;

% now feed Ad,Bd,Cd into the MPC builder from problem 8.1
N = 10;
Q = eye(size(Cd,1));
R = 0.1*eye(size(Bd,2));
mpc = mpc_unconstrained_design(Ad, Bd, Cd, N, Q, R);

xk = xs;  % start at steady state

for k = 1:200
    % compute input
    u_k = mpc_compute_unconstrained(mpc, xk, xs, u_s);

    % simulate one step of the discrete model
    xk = Ad * xk + Bd * u_k;   % (no disturbance here)
end


%% --- SIMULATE CLOSED-LOOP MPC AND PLOT ---

Tf = 200;                 % number of simulation steps
nx = size(Ad,1);
ny = size(Cd,1);
nu = size(Bd,2);

x = xs;                   % start at steady state
X = zeros(nx, Tf+1);
Y = zeros(ny, Tf+1);
U = zeros(nu, Tf);

X(:,1) = x;
Y(:,1) = Cd * x;

for k = 1:Tf
    % compute MPC input (unconstrained, from 8.3)
    u = mpc_compute_unconstrained(mpc, x, xs, u_s);

    % apply to discrete model
    x = Ad * x + Bd * u;

    % store
    X(:,k+1) = x;
    Y(:,k+1) = Cd * x;
    U(:,k)   = u;
end

t = 0:Tf;

mkdir('figures/problem_8');
% --------- PLOT HEIGHTS (outputs) ----------
figure;
%set figure size
set(gcf, 'Position', [100, 100, 800, 800]);
for i = 1:ny
    subplot(ny,1,i);
    plot(t, Y(i,:), 'LineWidth', 1.5);
    hold on;
    yref = Cd(i,:) * xs;           % steady-state height for output i
    yline(yref, '--k');
    grid on;
    xlabel('time k');
    ylabel(sprintf('h_%d [cm]', i));
    title(sprintf('Tank %d height', i));
end
sgtitle('MPC closed-loop heights');
% create new folder in figures folder for problem 8 to save plot as pdf
saveas(gcf, 'figures/problem_8/MPC_heights.pdf');

% --------- PLOT INPUTS ----------
figure;
for j = 1:nu
    subplot(nu,1,j);
    stairs(0:Tf-1, U(j,:), 'LineWidth', 1.5);
    hold on;
    uref = u_s(j);
    yline(uref, '--k');
    grid on;
    xlabel('time k');
    ylabel(sprintf('F_%d [cm^3/s]', j));
    title(sprintf('Pump %d flow', j));
end
sgtitle('MPC inputs');
saveas(gcf, 'figures/problem_8/MPC_inputs.pdf');