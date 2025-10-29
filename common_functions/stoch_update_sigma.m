%   Derivation for the stochastic added variables F3 and F4
function [d_diff,dw_it] = stoch_update_sigma(d,d_mean,dt,dw,dw_it)
    tau_v = [0.2; 0.1];
    sigma_d = [0.1; 0.1];

    d_diff = tau_v.*(d_mean - d)*dt + sigma_d.*(dw(:,dw_it));
    dw_it = dw_it + 1;
end