function [x,info] = Problem(H,g,l,u,A,bl,bu,xinit)
    %QPSOLVER for solving convex quadratic
        %
        %   min  phi = 1/2 * x'*H*x + g'*x
        %   s.t. l <= x <= u
        %        bl <= A*x <= b
        %   Inputs:
            %       H     : Hessian matrix (n x n, symmetric positive semi-definite)
            %       g     : Gradient vector (n x 1)
            %       l     : Lower bounds for x (n x 1 or scalar)
            %       u     : Upper bounds for x (n x 1 or scalar)
            %       A     : Linear inequality matrix (m x n)
            %       bl    : Lower bounds for A*x (m x 1)
            %       bu    : Upper bounds for A*x (m x 1)
            %       xinit : Initial guess for x (n x 1, optional)

    % Symmetric H
    H = (H + H') / 2;

        % Convert to quadprog's inequality form Aeq*x = beq, A*x <= b
            % bl <= Ax <= bu can be rewritten as:
            % Ax <= bu  (A_ineq = A, b_ineq = bu)
            % -Ax <= -bl (A_ineq = -A, b_ineq = -bl)
    
    % Combine into a single inequality system for quadprog
    A_ineq = [A; -A];
    b_ineq = [bu; -bl];

    % Set options for quadprog - verify with documentations if needed
    options = optimoptions('quadprog', ...
                           'Display', 'off', ... 
                           'Algorithm', 'interior-point-convex');

    % Call quadprog
    if nargin == 8 && ~isempty(xinit)
        [x, fval, exitflag, output, lambda] = quadprog(H, g, A_ineq, b_ineq, [], [], l, u, xinit, options);
    else
        [x, fval, exitflag, output, lambda] = quadprog(H, g, A_ineq, b_ineq, [], [], l, u, [], options);
    end

    % Populate info structure
    info.exitflag = exitflag;
    info.fval = fval;
    info.output = output;
    info.lambda = lambda;

    % Handler in case no solution was found
    if isempty(x)
        warning('qpsolver:NoSolution', 'quadprog did not find a solution.');
        x = NaN(size(g)); % Return NaN for x if no solution
    end

end