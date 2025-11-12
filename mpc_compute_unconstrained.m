function u = mpc_compute_unconstrained(mpc, xk, xs, u_s)
% mpc_compute_unconstrained
%   Compute MPC control move for current state xk
%   using the unconstrained MPC design in mpc (from 8.1).
%
%   mpc.Kmpc : (nu x nx)
%
%   xk  : current plant state (same as used in linearization)
%   xs  : steady-state state (same as used in fsolve)
%   u_s : steady-state input  (same as used in fsolve)
%
%   returns
%     u : actual input to send to the plant

    % deviation state
    x_tilde = xk - xs;

    % MPC move in deviation form
    u_tilde = mpc.Kmpc * x_tilde;

    % back to real input
    u = u_s + u_tilde;
end
