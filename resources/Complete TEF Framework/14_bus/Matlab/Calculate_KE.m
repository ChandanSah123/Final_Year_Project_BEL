function [KE_total, KE_corr] = Calculate_KE(npts, g, H, Ws, w, MOD)
    % Initialize output vectors
    KE_total = zeros(npts, 1);
    KE_corr  = zeros(npts, 1);
    
    % Temporary variable for relative speed calculation
    weq = zeros(npts, 1);

    for t = 1:npts
        
        % --- 1. Calculate Total System Kinetic Energy ---
        % Sum of K.E. of all individual generators
        for i = 1:g
            % Formula: 0.5 * M * w^2  where M = 2H/Ws
            KE_total(t) = KE_total(t) + 0.5 * (2 * H(i) / Ws) * (w(t,i))^2;
        end

        % --- 2. Calculate Corrected Kinetic Energy (MOD Based) ---
        wcr = 0;
        wsys = 0;
        Hcr = 0;
        Hsys = 0;
        
        for i = 1:g
            % Check if generator 'i' is in the MOD group
            is_critical = 0;
            for j = 1:numel(MOD)
                if i == MOD(j)
                    is_critical = 1;
                    break;
                end
            end
            
            % Sort into Critical Cluster (MOD) or Rest of System
            if is_critical == 1
                wcr = wcr + H(i) * w(t,i);
                Hcr = Hcr + H(i);
            else
                wsys = wsys + H(i) * w(t,i);
                Hsys = Hsys + H(i);
            end
        end
        
        % Calculate equivalent inertial center speeds
        if Hcr > 0, wcr = wcr / Hcr; end
        if Hsys > 0, wsys = wsys / Hsys; end
        
        % Calculate Equivalent Speed and Inertia
        weq(t) = wcr - wsys;
        
        if (Hcr + Hsys) > 0
            Heq = (Hcr * Hsys) / (Hcr + Hsys);
        else
            Heq = 0; % Avoid division by zero if inertias are invalid
        end
        
        % Store Corrected Kinetic Energy
        KE_corr(t) = 0.5 * (2 * Heq / Ws) * (weq(t))^2;
    end
end