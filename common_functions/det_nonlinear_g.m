%   Measurement function g(X,p)
function h = det_nonlinear_g(X,p)
    % Height h_i = m_i / (rho * A_i)
    h = [0;0;0;0];
    h(1) = X(1) / (p(12) * p(5));
    h(2) = X(2) / (p(12) * p(6));
    h(3) = X(3) / (p(12) * p(7));
    h(4) = X(4) / (p(12) * p(8));
    
end