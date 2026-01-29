function theta_u = Calculate_theta_u(MOD_idx, MOD_sort_data, num_gen, ths, H)
   ths = ths(:); 
    H   = H(:);
% Calculate M (Inertia constant)
    Ws = 2 * pi * 60;
    M = 2 * H / Ws;

    % Extract MOD Group
    MOD = MOD_sort_data(1:MOD_idx, 1);
    group1 = MOD;
    
    % Identify Group 2 (Rest of System)
    all_gen = 1:num_gen;
    group2 = setdiff(all_gen, group1);
    
    % Calculate Equivalent Inertias
    M1 = sum(M(group1));
    M2 = sum(M(group2));
    MT = M1 + M2;
    
    % Calculate COI of SEP angles
    theta1_s = sum(ths(group1) .* M(group1)) / M1;
    theta2_s = sum(ths(group2) .* M(group2)) / M2;
    
    % Calculate Angle Shift
    d_theta_s = theta1_s - theta2_s;
    d_theta1 = (pi - (2 * d_theta_s)) * (M2 / MT);
    d_theta2 = (pi - (2 * d_theta_s)) * (M1 / MT);
    
    % Initialize and Assign Output
    theta_u = zeros(num_gen, 1);
    theta_u(group1) = ths(group1) + d_theta1;     % belong to group1
    theta_u(group2) = ths(group2) - d_theta2;     % belong to group2
end