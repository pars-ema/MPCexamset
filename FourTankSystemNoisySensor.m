function y = FourTankSystemNoisySensor(x_k, p, Rvv)
% FOURTANKSYSTEMSENSOR  Returns noisy height measurements
%   y = FourTankSystemSensor(x_k, p, Rvv)
%   x_k  : state vector (4x1) or time-series matrix (N x 4)
%   p    : parameter vector (same as your p)
%   Rvv  : 4x1 vector of noise variances for each tank
%
%   y    : measurements (same size as x_k) in units of heights [cm]

% Extract tank areas and density
A = p(5:8);       % [cm^2]
rho = p(12);

% Convert masses -> heights
h = x_k ./ (rho * A');

% If Rvv not given, default to small variance
if nargin < 3 || isempty(Rvv)
    Rvv = 0.01*ones(1,4);   % default variances
end

% Generate Gaussian noise for each measurement
sigma = sqrt(Rvv);          % standard deviations
v = randn(size(h)) .* sigma; 

% Add noise
y = h + v;
size(y)
end