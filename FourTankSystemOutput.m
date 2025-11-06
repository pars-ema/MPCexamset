function y = FourTankSystemOutput(x_k,p)
%FOURTANKSYSTEMSENSOR Summary of this function goes here
%   Detailed explanation goes here

A = p(5:8,1)';
rho = p(12)
y = x_k./(rho*A);
y = y(1:2);

end

