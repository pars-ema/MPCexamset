%   The heights (h_x) are the outputs
function z = det_nonlinear_h(X,p)
    z = det_nonlinear_g(X,p);
    z = z(1:2);
end