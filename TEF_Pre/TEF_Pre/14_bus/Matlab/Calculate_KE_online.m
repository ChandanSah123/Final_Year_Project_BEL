function KE_corr = Calculate_KE_online(g, H, Ws, w, MOD)
    % Initialize output vector
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
                wcr = wcr + H(i) * w(i);
                Hcr = Hcr + H(i);
            else
                wsys = wsys + H(i) * w(i);
                Hsys = Hsys + H(i);
            end
        end
        
        % Calculate equivalent inertial center speeds
        if Hcr > 0, wcr = wcr / Hcr; end
        if Hsys > 0, wsys = wsys / Hsys; end
        
        % Calculate Equivalent Speed and Inertia
        weq = wcr - wsys;
        
        if (Hcr + Hsys) > 0
            Heq = (Hcr * Hsys) / (Hcr + Hsys);
        else
            Heq = 0; % Avoid division by zero if inertias are invalid
        end
        
        % Store Corrected Kinetic Energy
        KE_corr = 0.5 * (2 * Heq / Ws) * (weq)^2;
  
end