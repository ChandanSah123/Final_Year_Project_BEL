function x = NR(f_func, J_func, x0, accuracy, max_itr)
% Newton–Raphson with COI constraint and damping
%
% Solves f(x)=0 where f is power-angle equations
% Fixes singular Jacobian due to rotational invariance

x = x0;

if nargin < 4; accuracy = 1e-8; end
if nargin < 5; max_itr  = 40;   end

% ---- USER PARAMETERS ----
damping = 0.6;     % <--- critical for convergence
M_gen = evalin('base','M_gen');   % inertia vector
M_T   = sum(M_gen);

for iter = 1:max_itr

    f = f_func(x);
    J = J_func(x);

    % -------------------------------
    % COI CONSTRAINT (KEY FIX)
    % -------------------------------
    % Replace first equation with:
    % sum(M_i * theta_i) = 0
    J(1,:) = M_gen' / M_T;
    f(1)   = 0;

    % Solve linear system
    dx = -J \ f;

    % Damped update
    x = x + damping * dx;

    % Convergence test
    if norm(dx) < accuracy && norm(f) < accuracy
        return;
    end
end

disp('NR did not converge!');
x = inf;
end
