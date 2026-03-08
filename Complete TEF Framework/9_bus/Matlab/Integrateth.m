function dy = Integrateth(~, y, Pi, C, D, H)
    % INTEGRATETH Generic gradient calculation for any number of machines
    % Solves for the gradient: f_i(theta) = P_acc_i - (M_i/M_tot) * P_COI
    
    g = length(y); % Detect system size (e.g., 10)
    
    % 1. Calculate Electrical Power (Pe) & Accelerating Power (P_acc)
    % P_acc = Pi - Pe
    P_acc = zeros(g, 1);
    
    for i = 1:g
        Pe_i = 0;
        for j = 1:g
            if i ~= j
                theta_diff = y(i) - y(j);
                % Pe = sum( C*sin(dij) + D*cos(dij) )
                Pe_i = Pe_i + C(i,j)*sin(theta_diff) + D(i,j)*cos(theta_diff);
            end
        end
        P_acc(i) = Pi(i) - Pe_i;
    end
    
    % 2. Calculate P_COI (Total System Accelerating Power)
    % Sum of all individual accelerating powers
    P_COI_total = sum(P_acc);
    
    % 3. Calculate Gradient in COI Frame
    % dy = P_acc_i - (H_i / H_sum) * P_COI_total
    H_sum = sum(H);
    dy = P_acc - (H / H_sum) * P_COI_total;
    
end