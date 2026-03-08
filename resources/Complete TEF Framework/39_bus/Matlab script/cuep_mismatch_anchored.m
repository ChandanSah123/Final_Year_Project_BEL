function F_aug = cuep_mismatch_anchored(theta, Pm, E, Y, M, M_tot)
    % 1. Standard Mismatch (N equations)
    V = E .* (cos(theta) + 1i * sin(theta));
    I = Y * V;
    Pe = real(V .* conj(I));
    P_acc = Pm - Pe;
    P_COI = sum(P_acc);
    
    % Eq 1..N: Power Mismatch in COI Frame
    F_mismatch = P_acc - (M ./ M_tot) * P_COI;
    
    % 2. Anchor Constraint (1 equation)
    % The weighted sum of angles must be zero (definition of COI)
    % We scale this by 10 to give it "weight" in the optimization
    F_anchor = (sum(M .* theta) / M_tot) * 10; 
    
    % 3. Combine into (N+1) vector
    F_aug = [F_mismatch; F_anchor];
end