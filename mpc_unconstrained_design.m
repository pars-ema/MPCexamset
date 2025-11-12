function mpc = mpc_unconstrained_design(A, B, C, N, Q, R)
% mpc_unconstrained_design  Build matrices for unconstrained MPC
%
%   System:
%       x(k+1) = A x(k) + B u(k)
%       z(k)   = C x(k)
%
%   Cost:
%       J = 1/2 sum_{j=1}^N z_{k+j|k}' Q z_{k+j|k} + ...
%           1/2 sum_{j=0}^{N-1} u_{k+j|k}' R u_{k+j|k}
%
%   Inputs:
%       A,B,C : discrete-time state space matrices
%       N     : prediction horizon
%       Q,R   : stage weights
%
%   Output (struct mpc):
%       mpc.Phi     -- output prediction matrix (N*ny x nx)
%       mpc.Gamma   -- output-input prediction matrix (N*ny x N*nu)
%       mpc.Qbar    -- blockdiag(Q,...,Q)
%       mpc.Rbar    -- blockdiag(R,...,R)
%       mpc.H       -- Hessian of QP
%       mpc.F       -- state coupling matrix
%       mpc.Kmpc    -- first-move MPC gain (u = Kmpc * x)
%
%   This solves problem 8.1 (design part): it prepares the
%   unconstrained MPC matrices based on a discrete-time model.

    % dimensions
    nx = size(A,1);
    nu = size(B,2);
    ny = size(C,1);

    % ---------- build Phi (N*ny x nx) ----------
    Phi = zeros(N*ny, nx);
    A_power = A;
    for i = 1:N
        Phi((i-1)*ny+1:i*ny, :) = C * A_power;
        A_power = A_power * A;
    end

    % ---------- build Gamma (N*ny x N*nu) ----------
    Gamma = zeros(N*ny, N*nu);
    for row = 1:N
        % CA^{row-1} B goes into the (row,1) input block, etc.
        for col = 1:row
            A_block = A^(row-col);
            Gamma( (row-1)*ny+1:row*ny, (col-1)*nu+1:col*nu ) = C * A_block * B;
        end
    end

    % ---------- build block-diagonal Qbar, Rbar ----------
    Qbar = kron(eye(N), Q);
    Rbar = kron(eye(N), R);

    % ---------- build QP matrices ----------
    H = Gamma' * Qbar * Gamma + Rbar;   % (N*nu x N*nu)
    F = Gamma' * Qbar * Phi;            % (N*nu x nx)

    % ---------- compute first-move gain ----------
    % solve H * G = F  -> G = H \ F
    G = H \ F;  % (N*nu x nx)
    % select first nu rows (current control move)
    E = [eye(nu), zeros(nu, (N-1)*nu)];
    Kmpc = - E * G;  % (nu x nx)

    % ---------- pack output ----------
    mpc.A     = A;
    mpc.B     = B;
    mpc.C     = C;
    mpc.N     = N;
    mpc.Phi   = Phi;
    mpc.Gamma = Gamma;
    mpc.Qbar  = Qbar;
    mpc.Rbar  = Rbar;
    mpc.H     = H;
    mpc.F     = F;
    mpc.Kmpc  = Kmpc;
end

