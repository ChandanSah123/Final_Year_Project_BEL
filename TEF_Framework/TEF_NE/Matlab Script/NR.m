function x = NR(f_func, J_func, x0, accuracy, max_itr)
% Newton–Raphson with conditional COI constraint
% Works for:
%   1) SEP solving (theta-only)
%   2) Time-domain implicit integration (theta + omega)

x = x0;

if nargin < 4; accuracy = 1e-8; end
if nargin < 5; max_itr  = 50;   end

% Parameters from base workspace
M_gen = evalin('base','M_gen');
M_T   = sum(M_gen);

damping = 0.6;   % stable value

for iter = 1:max_itr

    f = f_func(x);
    J = J_func(x);

    % --------------------------------------------------
    % APPLY COI CONSTRAINT ONLY FOR THETA-ONLY PROBLEM
    % --------------------------------------------------
    if length(x) == length(M_gen)
        % Replace first equation with COI constraint
        J(1,:) = M_gen' / M_T;
        f(1)   = 0;
    end

    dx = -J \ f;

    x = x + damping * dx;

    if norm(dx) < accuracy && norm(f) < accuracy
        return;
    end
end

disp('NR did not converge!');
x = inf;
end
