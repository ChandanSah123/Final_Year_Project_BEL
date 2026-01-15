
function [dtheta1, dtheta2] =MOD_SHIFT(MOD_List, num_gen, thetas, H)
% GET_MOD_SHIFTS Calculates angle shifts for CUEP using MOD technique.
%
% Inputs:
%   MOD_List : Vector of generator IDs in the Critical Cluster (e.g., [2])
%   num_gen  : Total number of generators (integer)
%   thetas   : Vector of Stable Equilibrium Point (SEP) angles [rad]
%   H        : Vector of Inertia Constants for all generators [MWs/MVA]
%
% Outputs:
%   dtheta1  : Shift angle for the Critical Group (Group I)
%   dtheta2  : Shift angle for the Rest Group (Group II)

    %% 1. Identify Groups
    % Group I: Machines in the MOD list
    Group1 = MOD_List;
    
    % Group II: All other machines (Rest)
    all_gens = 1:num_gen;
    Group2 = setdiff(all_gens, Group1);
    
    %% 2. Calculate Group Inertias
    % M_I and M_II as per Equation (3)
    M1 = sum(H(Group1)); 
    M2 = sum(H(Group2)); 
    MT = M1 + M2;        % Total Inertia
    
    %% 3. Calculate Center of Inertia (COI) Angles
    % theta_I^s and theta_II^s as per Equation (3)
    
    % Weighted sum for Group I
    theta1_s = sum(thetas(Group1) .* H(Group1)) / M1;
    
    % Weighted sum for Group II
    theta2_s = sum(thetas(Group2) .* H(Group2)) / M2;
    
    %% 4. Calculate Angle Shifts
    % Separation between the two groups
    theta_sep = theta1_s - theta2_s; % (theta_I - theta_II)
    
    % Calculate shifts using Equation (4) from the paper
    % dtheta1 is for the Critical Group (Advanced)
    dtheta1 = (pi - 2 * theta_sep) * (M2 / MT);
    
    % dtheta2 is for the Rest Group
    dtheta2 = (pi - 2 * theta_sep) * (M1 / MT);

end