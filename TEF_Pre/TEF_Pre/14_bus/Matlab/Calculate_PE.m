function PE = Calculate_PE(npts, g, Pi, C, D, th, ths)
    % Initialize vectors
    PE1 = zeros(npts, 1);
    PE2 = zeros(npts, 1);
    PE3 = zeros(npts, 1);
    PE  = zeros(npts, 1);
    
    for t = 1:npts
        
        % Term 1: Position Energy
        for i = 1:g
            % FIXED: Added NEGATIVE sign
            % (Standard formula: - Sum Pi * (theta - theta_s))
            PE1(t) = PE1(t) - Pi(i) * (th(t,i) - ths(i));
        end
        
        % Terms 2 & 3
        for i = 1:g-1
            for j = i+1:g
                th_ij = th(t,i) - th(t,j);
                ths_ij = ths(i) - ths(j);
                
                % Term 2: Magnetic Energy
                % FIXED: Added NEGATIVE sign
                % (Standard formula: - Sum C * (cos(th) - cos(ths)))
                % Since (cos_th - cos_ths) is usually negative as you move, 
                % this double negative makes the Energy POSITIVE.
                PE2(t) = PE2(t) - C(i,j) * (cos(th_ij) - cos(ths_ij));
                
                % Term 3: Dissipative Energy
                % (This sign depends on your D definition, but usually opposes motion)
                term_num = (th(t,i) - ths(i) + th(t,j) - ths(j));
                term_den = (th_ij - ths_ij);
                
                if abs(term_den) < 1e-9
                    ratio = 0;
                else
                    ratio = term_num / term_den;
                end
                
                % Note: If you flip the others, verify if D needs flipping. 
                % Usually, D is path dependent addition.
                PE3(t) = PE3(t) + D(i,j) * ratio * (sin(th_ij) - sin(ths_ij));
            end
        end
        
        % Total Potential Energy
        PE(t) = PE1(t) + PE2(t) + PE3(t);
    end
end