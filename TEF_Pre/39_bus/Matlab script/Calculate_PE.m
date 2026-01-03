function PE = Calculate_PE(npts, g, Pi, C, D, th, ths)
    % Calculates Standard Potential Energy for the whole trajectory
    PE1 = zeros(npts, 1);
    PE2 = zeros(npts, 1);
    PE3 = zeros(npts, 1);
    PE  = zeros(npts, 1);
    
    for t = 1:npts
        % Term 1: Position Energy ( - sum Pi * theta_diff )
        for i = 1:g
            PE1(t) = PE1(t) - Pi(i) * (th(t,i) - ths(i));
        end
        
        % Terms 2 & 3
        for i = 1:g-1
            for j = i+1:g
                th_ij  = th(t,i) - th(t,j);
                ths_ij = ths(i) - ths(j);
                
                % Term 2: Magnetic ( - C * (cos - cos) )
                PE2(t) = PE2(t) - C(i,j) * (cos(th_ij) - cos(ths_ij));
                
                % Term 3: Dissipative ( + Integral D )
                term_num = (th(t,i) - ths(i) + th(t,j) - ths(j));
                term_den = th_ij - ths_ij;
                
                if abs(term_den) < 1e-9
                    ratio = 0;
                else
                    ratio = term_num / term_den;
                end
                
                PE3(t) = PE3(t) + D(i,j) * ratio * (sin(th_ij) - sin(ths_ij));
            end
        end
        
        PE(t) = PE1(t) + PE2(t) + PE3(t);
    end
end