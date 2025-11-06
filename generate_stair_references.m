function [r1, r2] = generate_stair_references(t, t_change_1, values_1, t_change_2, values_2)
%GENERATE_STAIR_REFERENCES Generate piecewise constant (stair) reference signals
%
% [r1, r2] = generate_stair_references(t, t_change_1, values_1, t_change_2, values_2)
%
% This function creates two piecewise-constant reference trajectories 
% (for example, desired tank levels) that follow user-defined step times 
% and values. Each reference will have the same length as the time vector `t`.
%
% INPUTS:
%   t            : time vector [1 x N] (e.g. 0:1:3600)
%   t_change_1   : vector of times [s] when r1 changes value
%   values_1     : vector of r1 values (length must match t_change_1)
%   t_change_2   : vector of times [s] when r2 changes value
%   values_2     : vector of r2 values (length must match t_change_2)
%
% OUTPUTS:
%   r1, r2       : vectors [1 x N] with piecewise-constant reference signals
%
% Example:
%   t = 0:1:3600;
%   [r1, r2] = generate_stair_references(t, [0 1200 2400], [10 15 12], ...
%                                        [0 1800], [8 9]);
%
% Author: MPC Project Team
% -------------------------------------------------------------

    % --- Safety checks ---
    if length(t_change_1) ~= length(values_1)
        error('t_change_1 and values_1 must have the same length');
    end
    if length(t_change_2) ~= length(values_2)
        error('t_change_2 and values_2 must have the same length');
    end

    % --- Initialize outputs ---
    r1 = zeros(size(t));
    r2 = zeros(size(t));

    % --- Generate first reference ---
    for i = 1:length(t_change_1)-1
        idx = (t >= t_change_1(i)) & (t < t_change_1(i+1));
        r1(idx) = values_1(i);
    end
    r1(t >= t_change_1(end)) = values_1(end);

    % --- Generate second reference ---
    for i = 1:length(t_change_2)-1
        idx = (t >= t_change_2(i)) & (t < t_change_2(i+1));
        r2(idx) = values_2(i);
    end
    r2(t >= t_change_2(end)) = values_2(end);
end
