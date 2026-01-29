function theta_u = Calculate_theta_u_online1( MOD, num_gen, ths, ~) 
  
    group1 = MOD;
    
    % Identify Group 2 (Rest of System)
    all_gen = 1:num_gen;
    group2 = setdiff(all_gen, group1);
     
    % Initialize and Assign Output
    theta_u = zeros(num_gen, 1);
    theta_u(group1)=pi-ths(group1);     % belong to group1
    theta_u(group2) =ths(group2);   % belong to group2
end