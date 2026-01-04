function f_total = Calculate_Fsum(th, g, Pi, C, D, H)
    % Calculates the sum of ABSOLUTE accelerating powers (fsum)
    % exactly matching your original logic.
    
    Pcoi = 0;
    % 1. Calculate COI Accelerating Power (Global System Drift)
    for i = 1:g
        Pcoi = Pcoi + Pi(i);
        for j = i+1:g
            % Summing electrical interactions for COI
            Pcoi = Pcoi - 2*D(i,j)*cos(th(i) - th(j)); 
        end
    end
    
    M_tot = sum(H); 
    f_total = 0;
    
    % 2. Calculate Individual Mismatches
    for i = 1:g
        Pei = 0;
        for j = 1:g
            if j ~= i
                Pei = Pei + C(i,j)*sin(th(i) - th(j)) + D(i,j)*cos(th(i) - th(j));
            end
        end
        
        % The mismatch equation in COI frame:
        f_val = Pi(i) - Pei - (H(i)/M_tot) * Pcoi;
        
        % CORRECTION: Sum of Absolute Values (L1 Norm)
        f_total = f_total + abs(f_val); 
    end
end