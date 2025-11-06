% This script illustrate creation of a discrete-time state space model
% from a continuous-time MIMO transfer function. It uses the function
% mimoctf2dss. 
clear
clc


%% Create transfer function matrix representation G(s)
%  see slides
num = cell(2,2);
den = cell(2,2);
lambda = zeros(2,2);

% G(1,1)
num(1,1) = {0.1133}
den(1,1) = {[1.783 4.48 1]}
lambda(1,1) = 0.715

% G(2,1)
num(2,1) = {0.3378}
den(2,1) = {[0.361 1.09 1]}
lambda(2,1) = 0.299

% G(1,2)
num(1,2) = {0.9222}
den(1,2) = {[2.071 1]}
lambda(1,2) = 0.0

% G(2,2)
num(2,2) = {-0.321}
den(2,2) = {[0.104 2.463 1]}
lambda(2,2) = 0.9

%% Create the discrete-time state space model
Ts   = 0.05     % sample time
Nmax = 100      % maximum state dimension is Nmax/2
tol  = 1.0e-8   % tolerance

[Ad,Bd,Cd,Dd,sH] = mimoctf2dss(num,den,lambda,Ts,Nmax,tol);

whos 
%%
figure(1)
semilogy(sH,'b.-','linewidth',2,'markersize',30)