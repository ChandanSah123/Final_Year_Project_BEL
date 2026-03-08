function [dy, Jac] = Integrateth(~, y, Pi, C, D, H)
    % INTEGRATETH_JAC Calculates Mismatch AND Jacobian
    g = length(y);
    
    % --- 1. Calculate Mismatch (dy) ---
    P_acc = zeros(g, 1);
    
    % We also need J_acc (Partial derivatives of P_acc)
    % J_acc(i, j) = d(P_acc_i) / d(theta_j)
    J_acc = zeros(g, g); 
    
    for i = 1:g
        Pe_i = 0;
        for j = 1:g
            if i ~= j
                theta_diff = y(i) - y(j);
                sin_t = sin(theta_diff);
                cos_t = cos(theta_diff);
                
                % Power Term
                Pe_i = Pe_i + C(i,j)*sin_t + D(i,j)*cos_t;
                
                % Jacobian Terms (Derivatives)
                % d(P_acc_i)/d(theta_j) = - d(Pe_i)/d(theta_j)
                % d(Pe_i)/d(theta_j) involves -1 from chain rule (theta_i - theta_j)
                % Result: C*cos - D*sin
                deriv_term = C(i,j)*cos_t - D(i,j)*sin_t;
                
                J_acc(i, j) = deriv_term; % Off-diagonal
            end
        end
        P_acc(i) = Pi(i) - Pe_i;
        
        % Diagonal element is negative sum of off-diagonals
        % (Because changing theta_i affects every pair involving i)
        J_acc(i, i) = -sum(J_acc(i, :)); 
    end
    
    % --- 2. COI Transformation for Mismatch ---
    P_COI_total = sum(P_acc);
    H_sum = sum(H);
    h_ratio = H / H_sum; % Vector (g x 1)
    
    dy = P_acc - h_ratio * P_COI_total;
    
    % --- 3. COI Transformation for Jacobian ---
    % The system function is F = P_acc - h_ratio * sum(P_acc)
    % The Jacobian is J = J_acc - h_ratio * sum(rows of J_acc)
    
    % sum(J_acc, 1) sums the columns, effectively d(P_COI)/d(theta_j)
    dP_COI_dtheta = sum(J_acc, 1); 
    
    % Broadcast subtraction
    Jac = J_acc - (h_ratio * dP_COI_dtheta);

end