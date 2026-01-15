function dy = Integrateth2(~, y)
    %% 1. Define Parameters (Fill these in)
    global Yint_post;
    global data;% Assuming this is available
    
    % PLACEHOLDERS (You add the values later)
    H = [23.64;6.4;3.01]; 
    E=[1.054 1.050 1.017];
    Pm = [0.716410160064697	1.62999987602234	0.850000023841858];
    Y1=Yint_post;
    
    
    %% 2. Pre-calculate Matrices C and D
    % The original code calculated C and D inside loops. We can do it in one line.
    % We assume E is a column vector (g x 1).
    
    g = length(E);
    
    % Create a matrix where Element (i,j) = E(i) * E(j)
    EE_matrix = E * E'; 
    
    % Calculate C (Imaginary/Sine component) and D (Real/Cosine component)
    % .* does element-wise multiplication
    C = EE_matrix .* imag(Y1);
    D = EE_matrix .* real(Y1);
    
    %% 3. Calculate Angle Differences Matrix
    % We need (theta_i - theta_j) for all combinations.
    % This creates a (g x g) matrix where element (i,j) is y(i) - y(j)
    Theta_diff = y - y'; 
    
    %% 4. Calculate Electrical Power (Pe) using C and D
    % The formula: Sum( C_ij * sin(theta_ij) + D_ij * cos(theta_ij) )
    % sum(..., 2) sums across the rows (j=1 to g) for each generator i
    
    Pe = sum( C .* sin(Theta_diff) + D .* cos(Theta_diff), 2 );
    
    %% 5. Calculate Accelerating Power (COI Frame)
    % This part remains the same logic as the original code
    
    Total_H = sum(H);
    P_acc   = Pm - Pe; % Accelerating power in local frame
    
    % COI Correction: P_acc_COI = P_acc - (H_i/H_T) * Sum(P_acc)
    dy = P_acc - (H / Total_H) * sum(P_acc);
    
end