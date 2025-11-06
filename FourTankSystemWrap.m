function xdot = FourTankSystemWrap(x, u, d, p)
%FOURTANKSYSTEMWRAP Summary of this function goes here
%   Detailed explanation goes here
xdot = FourTankSystemModified(0,x,u,d,p);
end

