function [theta_sol, fval] = Solve_Equilibrium(x0,g,Pi,C,D ,H, type)
    % Common setup
    num_gen = length(Pi);
    
    % Method A approach (Recommended for both)
    % We define the objective as minimizing the power mismatch
    obj_fun = @(x)  Calculate_Fsum(x, g, Pi, C, D, H); % Assuming this function exists
    
    % COI Constraints (Always good to enforce for stability studies)
    Aeq = H';       
    beq = 0;
    
    % Bounds (Relaxed bounds for general search)
    lb = -2*pi * ones(num_gen, 1);
    ub =  2*pi * ones(num_gen, 1);
    
    opts = optimset('Algorithm', 'interior-point', 'Display', 'off');

    % Solve
    fprintf('Solving for %s...\n', type);
    [x_sol, fval] = fmincon(obj_fun, x0, [], [], Aeq, beq, lb, ub, [], opts);
    
    % Post-process: Center the angles (COI Correction)
    theta_sol = x_sol - (sum(x_sol .* H(:)) / sum(H(:)));
    
    % Sanity Check
    if fval > 1e-4
        warning('Solver did not converge to a true equilibrium! Residual error: %f', fval);
    else
        fprintf('%s found successfully.\n', type);
    end
end