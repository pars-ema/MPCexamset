%   The sigma noise is given as a diagonal matrix/format
function [dX,d_diff,dw_it] = sde_nonlinear_f(X,u,p,d_state,d_mean,dt,dw,dw_it)
    dX = stoch_nonlinear_f(X,u,d_state,p);
    [d_diff,dw_it] = stoch_update_sigma(d_state,d_mean,dt,dw,dw_it);
end