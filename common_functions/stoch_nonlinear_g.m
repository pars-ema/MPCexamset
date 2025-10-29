%   2. Measurement function g(x,p)
%   The measurements are the liquid heights (h_i) in cm.
%   y(t) = g(x(t),p) + v(t)
function h = stoch_nonlinear_g(X,p,R)
    %   For this application, R should be a vector 1x4, representing a
    %   diagonal matrix.
    v = [normrnd(0, R(1));normrnd(0, R(2));normrnd(0, R(3));normrnd(0, R(4))];

    % Height h_i = m_i / (rho * A_i)
    h = [0;0;0;0];
    h(1) = (X(1) / (p(12) * p(5))) + v(1);
    h(2) = (X(2) / (p(12) * p(6))) + v(2);
    h(3) = (X(3) / (p(12) * p(7))) + v(3);
    h(4) = (X(4) / (p(12) * p(8))) + v(4);
    
end