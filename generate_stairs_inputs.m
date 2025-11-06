function [u1, u2] = generate_stair_inputs(t, t_change_1, values_1, t_change_2, values_2)
%GENERATE_STAIR_INPUTS Generate piecewise constant (stair) input signals
%
% [u1, u2] = generate_stair_inputs(t, t_change_1, values_1, t_change_2, values_2)
%
% This function creates two piecewise-constant input vectors u1(t) and u2(t)
% that follow user-defined step times and values. Each input will have the
% same length as the time vector `t`.
%
% INPUTS:
%   t            : time vector [1 x N] (e.g. 0:1:3600)
%   t_change_1   : vector of times [s] when u1 changes value
%   values_1     : vector of u1 values (length must be length(t_change_1))
%   t_change_2   : vector of times [s] when u2 changes value
%   values_2     : vector of u2 values (length must be length(t_change_2))
%
% OUTPUTS:
%   u1, u2       : vectors [1 x N] with piecewise-constant signals

% --- Safety checks ---
if length(t_change_1) ~= length(values_1)
    error('t_change_1 and values_1 must have the same length');
end
if length(t_change_2) ~= length(values_2)
    error('t_change_2 and values_2 must have the same length');
end

% --- Initialize outputs ---
u1 = zeros(size(t));
u2 = zeros(size(t));

% --- Generate first input ---
for i = 1:length(t_change_1)-1
    idx = (t >= t_change_1(i)) & (t < t_change_1(i+1));
    u1(idx) = values_1(i);
end
% Last segment
u1(t >= t_change_1(end)) = values_1(end);

% --- Generate second input ---
for i = 1:length(t_change_2)-1
    idx = (t >= t_change_2(i)) & (t < t_change_2(i+1));
    u2(idx) = values_2(i);
end
u2(t >= t_change_2(end)) = values_2(end);

end
