function f_total = Calculate_Fsum(th, g, Pi, C, D, H)
    % 1. Compute Electrical Power Output (Pe) for ALL generators
    Pe = zeros(g, 1);
    for i = 1:g
        for j = 1:g
            if i ~= j
                Pe(i) = Pe(i) + C(i,j)*sin(th(i)-th(j)) + D(i,j)*cos(th(i)-th(j));
            end
        end
    end
    
    % 2. Calculate Raw Accelerating Power (P_acc)
    % Force Pi to be a column vector to prevent accidental matrix expansion
    Pi_col = Pi(:); 
    P_acc = Pi_col - Pe; 
    
    % 3. Calculate COI Accelerating Power
    P_COI_Total = sum(P_acc);
    M_tot = sum(H);
    
    % 4. Calculate Mismatch in COI Frame
    f_total = 0;
    for i = 1:g
        mismatch = P_acc(i) - (H(i) / M_tot) * P_COI_Total;
        f_total = f_total + mismatch.^2; % Sum of squares objective
    end
end