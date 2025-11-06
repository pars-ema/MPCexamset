function xdot = FourTankSystemModified(t, x, u, d, p)
% FOURTANKSYSTEM Model dx/dt = f(t,x,u,p) for 4-tank system
%
% This function implements a differential equation model
% for the 4-tank system.
%
% Syntax: xdot = FourTankSystem(t,x,u,p)

% Unpack states, MVs, and parameters
m = x;                         % Mass of liquid in each tank [g]
F = u;                         % Flow rates in pumps [cm3/s]
a = p(1:4);                    % Pipe cross sectional areas [cm2]
A = p(5:8);                    % Tank cross sectional areas [cm2]
g = p(9);                      % Acceleration of gravity [cm/s2]
gamma = p(10:11);              % Valve positions [-]
rho = p(12);                   % Density of water [g/cm3]
F3 = d(1);                      % Flow rate as disturbance to tank 3 [g/cm3]
F4 = d(2);                      % Flow rate as disturbance to tank 4 [g/cm3]


% Inflows
qin = zeros(4,1);
qin(1) = gamma(1) * F(1);      % Valve 1 to tank 1 [cm3/s]
qin(2) = gamma(2) * F(2);      % Valve 2 to tank 2 [cm3/s]
qin(3) = (1 - gamma(2)) * F(2) + F3;% Valve 2 to tank 3 [cm3/s]
qin(4) = (1 - gamma(1)) * F(1) + F4;% Valve 1 to tank 4 [cm3/s]

% Outflows
h = m ./ (rho * A);            % Liquid level in each tank [cm]
qout = a .* sqrt(max(0, 2 * g * h)); % Prevent negative heights

% Differential equations, mass balances
xdot = zeros(4,1);
xdot(1) = rho * (qin(1) + qout(3) - qout(1)); % Tank 1
xdot(2) = rho * (qin(2) + qout(4) - qout(2)); % Tank 2
xdot(3) = rho * (qin(3) - qout(3) + F3);           % Tank 3
xdot(4) = rho * (qin(4) - qout(4) + F4);           % Tank 4
