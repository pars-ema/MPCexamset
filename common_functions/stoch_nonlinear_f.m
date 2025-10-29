%   (F3,F4) are now defined and included in the scope of this exercise.
%   Those should be defined stochastic but piecewise constant

%   d = [F3; F4] (stochastic flow disturbances in cm3/s)
function dX = stoch_nonlinear_f(X,u,d,p)
    X_c = max(0, X);                    %   Avoid negative square roots
    dX = det_nonlinear_f(X_c,u,p);
    
    dX(3) = dX(3) + p(12)*d(1);
    dX(4) = dX(4) + p(12)*d(2);
end