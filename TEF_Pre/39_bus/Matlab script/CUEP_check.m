%% CHECK EQUILIBRIUM VALIDITY
% 1. Define the MOD angle vector you want to test
% assigning the result from your minimization to theta_MOD for the check
theta_MOD = theta_cuep_mod; 
% Ensure it is a column vector for consistency
if size(theta_MOD, 2) > 1; theta_MOD = theta_MOD'; end

% 2. Ensure the global admittance matrix is accessible
global Yint_post;
if isempty(Yint_post)
    % Not critical for the calculation below (which uses C/D/Pi), 
    % but kept for consistency with your template.
    warning('Global variable Yint_post is empty. Make sure your environment is set up.');
end

% 3. Call the function to calculate accelerating power (Gradient)
% CORRECTION: Passed Pi, C, D, H because Integrateth definition requires them.
f_val = Integrateth(0, theta_MOD, Pi, C, D, H);

% 4. Display and Analyze Results
disp('========================================');
disp('      EQUILIBRIUM POINT VERIFICATION    ');
disp('========================================');

fprintf('Testing Angles (rad) [First 3]: [%.4f; %.4f; %.4f]...\n', theta_MOD(1:3));

disp('----------------------------------------');
disp('Calculated Accelerating Power f(theta) (Mismatch):');
disp(f_val);

% Calculate the magnitude of the gradient vector
gradient_norm = norm(f_val);
fprintf('Gradient Norm: %.6f\n', gradient_norm);
disp('----------------------------------------');

% 5. Conclusion Logic
if gradient_norm < 1e-3
    disp('CONCLUSION: This IS a valid Equilibrium Point.');
    disp('The solver converged correctly (Gradient ~ 0).');
else
    disp('CONCLUSION: This is NOT a valid Equilibrium Point.');
    disp('The solver stopped prematurely or failed to converge (Gradient is High).');
end