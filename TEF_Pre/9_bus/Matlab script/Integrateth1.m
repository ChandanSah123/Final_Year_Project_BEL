function dy = Integrateth1(t, y)
    %% 1. Retrieve System Parameters
    % Define these variables or load them from your workspace/globals
    % global Y1;  % Un-comment if using global
    
    % PLACEHOLDERS: Add your specific data here
    % H  = [ ... ];  % Inertia Constants
    % E  = [ ... ];  % Internal Voltages
    % Pm = [ ... ];  % Mechanical Power
    % Y1 = [ ... ];  % Reduced Admittance Matrix

    %% 2. Prepare Vectors
    % Ensure all vectors are column vectors (N x 1) for matrix math
    y  = y(:);   % Generator Angles (theta)
    E  = E(:);   % Voltage Magnitudes
    H  = H(:);   % Inertia
    Pm = Pm(:);  % Mechanical Power

    %% 3. Vectorized Electrical Power Calculation
    % Instead of looping i and j, we use matrix operations.
    % V = E * e^(j * theta)
    V = E .* exp(1j * y); 
    
    % Current Injections: I = Y_bus * V
    I = Y1 * V;
    
    % Electrical Power: Pe = Real(V * conj(I))
    Pe = real(V .* conj(I));

    %% 4. Calculate Accelerating Power (Local Frame)
    % P_acc = Mechanical Power - Electrical Power
    P_acc = Pm - Pe;

    %% 5. Transform to COI Frame (Center of Inertia)
    % Formula: f_i = P_acc_i - (H_i / H_total) * sum(P_acc_all)
    
    Total_H     = sum(H);
    Total_P_acc = sum(P_acc);
    
    % Subtract the system-wide average acceleration weighted by inertia
    dy = P_acc - (H / Total_H) * Total_P_acc;

end