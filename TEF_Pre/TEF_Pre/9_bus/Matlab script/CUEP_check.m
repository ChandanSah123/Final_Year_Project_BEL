%% CHECK EQUILIBRIUM VALIDITY
% 1. Define the MOD angle vector you want to test
theta_MOD = [ -0.0255    0.2310    0.2918    0.3037    0.4291    0.3618    0.5286    0.2540 0.6277   -0.2603];
 
% theta_MOD= [0.5162 1.2452 0.9803];

% 2. Ensure the global admittance matrix is accessible
global Yint_post;
if isempty(Yint_post)
    error('Global variable Yint_post is empty. Run your main simulation first.');
end

% 3. Call the function to calculate accelerating power (Gradient)
% We pass 0 for the time argument (~) as it is ignored inside the function.
f_val = Integrateth(0, theta_MOD);

% 4. Display and Analyze Results
disp('========================================');
disp('      EQUILIBRIUM POINT VERIFICATION    ');
disp('========================================');
fprintf('Testing Angles (rad): [%.4f; %.4f; %.4f]\n', theta_MOD);
disp('----------------------------------------');
disp('Calculated Accelerating Power f(theta):');
disp(f_val);

% Calculate the magnitude of the gradient vector
gradient_norm = norm(f_val);
fprintf('Gradient Norm: %.6f\n', gradient_norm);
disp('----------------------------------------');

% 5. Conclusion Logic
if gradient_norm < 1e-3
    disp('CONCLUSION: This IS a valid Equilibrium Point.');
    disp('The solver converged correctly, but to a different/wrong saddle point.');
else
    disp('CONCLUSION: This is NOT a valid Equilibrium Point.');
    disp('The solver stopped prematurely or failed to converge.');
end