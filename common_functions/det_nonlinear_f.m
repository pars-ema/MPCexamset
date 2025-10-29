%   Prediction function f(X,u,p)
%   x = [m1; m2; m3; m4] (masses in g)
%   u = [F1; F2] (pump flow rates in cm3/s)
function dX = det_nonlinear_f(X,u,p)
    dX = [0;0;0;0];
    dX(1) =  p(12)*(p(9)*u(1)) + p(12)*p(3)*sqrt(2*p(11)*(X(3)/(p(12)*p(7)))) ...
            - p(12)*p(1)*sqrt(2*p(11)*(X(1)/(p(12)*p(5))));
    dX(2) =  p(12)*(p(10)*u(2)) + p(12)*p(4)*sqrt(2*p(11)*(X(4)/(p(12)*p(8)))) ...
            - p(12)*p(2)*sqrt(2*p(11)*(X(2)/(p(12)*p(6))));
    dX(3) =  p(12)*((1-p(10))*u(2)) - p(12)*p(3)*sqrt(2*p(11)*(X(3)/(p(12)*p(7))));
    dX(4) =  p(12)*((1-p(9))*u(1)) - p(12)*p(4)*sqrt(2*p(11)*(X(4)/(p(12)*p(8))));
end