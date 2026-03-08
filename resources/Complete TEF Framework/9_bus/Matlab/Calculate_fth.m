function fth = Calculate_fth(ths, th, Pi, npts, g, C, D, H)
    % Initialize output vector
    fth = zeros(npts, 1);
    
    % Initialize temporary vectors
    Pei = zeros(1, g);
    f   = zeros(1, g);
    
    for t = 1:npts
        % 1. Calculate Pcoi (Power of Center of Inertia)
        Pcoi = 0;
        for i = 1:g
            Pcoi = Pcoi + Pi(i);
            for j = i+1:g
                Pcoi = Pcoi - 2 * D(i,j) * cos(th(t,i) - th(t,j));
            end
        end
        % 2. Calculate Accelerating Power f(i) and fth dot product
        for i = 1:g
            Pei(i) = 0;
            for j = 1:g
                if j ~= i
                    Pei(i) = Pei(i) + C(i,j) * sin(th(t,i) - th(t,j)) + D(i,j) * cos(th(t,i) - th(t,j));
                end
            end
            
            % Accelerating power in COI frame
            f(i) = Pi(i) - Pei(i) - (H(i) / sum(H)) * Pcoi;
            
            % Accumulate fth (dot product of angle deviation and accel power)
            fth(t) = fth(t) + ((th(t,i) - ths(i)) * f(i));
        end
    end
end