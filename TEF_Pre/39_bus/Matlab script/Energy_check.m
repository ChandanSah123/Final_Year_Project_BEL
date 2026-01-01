% 1. Define your two results
%th_bcu = [-0.0255; 0.2310; 0.2918; 0.3037; 0.4291; 0.3618; 0.5286; 0.2540; 0.6277; -0.2603];
%th_mod = [-0.5700; -0.3135; -0.2527; -0.2409; -0.1154; -0.1828; -0.0159; -0.2905; 6.3664; -0.8048];
th_bcu=theta_cuep_pebs;
th_mod=theta_cuep_mod;

% 2. Normalize MOD (Wrap angles to -pi to +pi) to see if they match
th_mod_wrapped = mod(th_mod + pi, 2*pi) - pi;

fprintf('Machine 9 Comparison:\n');
fprintf('  BCU Angle: %.4f rad\n', th_bcu(9));
fprintf('  MOD Angle: %.4f rad (Raw)\n', th_mod(9));
fprintf('  MOD Angle: %.4f rad (Wrapped)\n', th_mod_wrapped(9));

% 3. Calculate Energy for both (Lowest Energy Barrier is Critical)
% Ensure Pi, C, D, g are in your workspace!
V_bcu = Calculate_PE_single_point(th_bcu, ths, Pi, C, D, g);
V_mod = Calculate_PE_single_point(th_mod, ths, Pi, C, D, g);

fprintf('\nPotential Energy (Energy Barrier):\n');
fprintf('  V_critical (BCU): %.4f\n', V_bcu);
fprintf('  V_critical (MOD): %.4f\n', V_mod);

if V_bcu < V_mod
    disp('RECOMMENDATION: Choose BCU. It represents the closer (more critical) instability barrier.');
else
    disp('RECOMMENDATION: Check Validity. MOD has lower energy, but might be numerically wrapped.');
end