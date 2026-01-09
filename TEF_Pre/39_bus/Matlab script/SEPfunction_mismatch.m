function F = SEPfunction_mismatch(x, Pm, E, C, D, H, Y)
    % 1. Robust Size Detection
    g = numel(x); 
    
    % 2. Calculate COI Term (P_coi)
    % This represents the total system loss correction
    Pcoi = 0;
    for i = 1:g 
        for j = i+1:g
            Pcoi = Pcoi + D(i,j) * cos(x(i)-x(j));
        end
    end
    Pcoi = Pcoi * 2;
    
    % 3. Calculate Mismatch Vector
    F = zeros(g, 1);
    
    % Pre-calculate constant term to save speed inside loop
    % (Weighted sum of losses)
    total_loss_term = sum(Pm - E.^2 .* real(diag(Y))');
    
    for i = 1:g
        % A. Standard Power Flow Mismatch
        for j = 1:g        
            if i ~= j
                F(i) = F(i) - (C(i,j)*sin(x(i)-x(j)) + D(i,j)*cos(x(i)-x(j)));
            else
                F(i) = F(i) + Pm(i) - E(i)^2 * real(Y(i,i));
            end
        end
        
        % B. Subtract COI Reference Frame Adjustment
        % (Distribute the system mismatch according to inertia H)
        F(i) = F(i) - (H(i)/sum(H)) * total_loss_term + (H(i)/sum(H)) * Pcoi;
    end
end