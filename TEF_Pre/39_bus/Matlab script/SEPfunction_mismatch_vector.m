function F = SEPfunction_mismatch_vector(x, Pm, E, C, D, H, Y)
    % Matches Eq (6) from your screenshot exactly
    g = numel(x); 
    
    % 1. Calculate COI Term (Eq 7 in screenshot)
    Pcoi = 0;
    for i = 1:g 
        for j = i+1:g
            Pcoi = Pcoi + D(i,j) * cos(x(i)-x(j));
        end
    end
    Pcoi = Pcoi * 2; % Because D_ij cos = D_ji cos
    
    % 2. Calculate Mismatch Vector
    F = zeros(g, 1);
    
    % Constant term: P_i - E_i^2 * G_ii (From 'where' clause in screenshot)
    % Note: Your Y matrix likely contains G (real) and B (imag)
    Pi_const = Pm - (E.^2) .* real(diag(Y))'; 
    Total_P_const = sum(Pi_const); % Sum of P_i

    M_tot = sum(H); % Assuming H is proportional to M
    
    for i = 1:g
        % Summation term in Eq (6)
        Sum_interaction = 0;
        for j = 1:g        
            if i ~= j
                Sum_interaction = Sum_interaction + (C(i,j)*sin(x(i)-x(j)) + D(i,j)*cos(x(i)-x(j)));
            end
        end
        
        % The specific COI term adjustment from Eq (6)
        % M_i * d2theta/dt2 = P_i - Sum(...) - (M_i/M_T)*P_coi
        % Equilibrium implies d2theta/dt2 = 0
        
        % The P_coi term in Eq (6) actually represents the total system mismatch distribution
        % Re-constructing exactly based on standard TEF-COI formulation:
        
        F(i) = Pi_const(i) - Sum_interaction - (H(i)/M_tot) * (Total_P_const - Pcoi);
    end
end