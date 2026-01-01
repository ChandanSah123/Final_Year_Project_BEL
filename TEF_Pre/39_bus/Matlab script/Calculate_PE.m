function PE =Calculate_PE(npts, g, Pi, C, D, th, ths)
    % Initialize vectors to store energy terms
    PE1 = zeros(npts, 1);
    PE2 = zeros(npts, 1);
    PE3 = zeros(npts, 1);
    PE  = zeros(npts, 1);

    for t = 1:npts
        % Term 1: Position Energy
        for i = 1:g
            PE1(t) = PE1(t) + Pi(i) * (th(t,i) - ths(i));
        end

        % Terms 2 & 3: Magnetic and Dissipative Energy
        for i = 1:g-1
            for j = i+1:g
                % Magnetic Energy
                PE2(t) = PE2(t) + C(i,j) * (cos(th(t,i) - th(t,j)) - cos(ths(i) - ths(j)));
                
                % Dissipative Energy (Linear Path)
                term_num = (th(t,i) - ths(i) + th(t,j) - ths(j));
                term_den = (th(t,i) - th(t,j) - ths(i) + ths(j));
                
                % Avoid division by zero if needed, otherwise use formula directly
                if abs(term_den) < 1e-9
                    ratio = 0;
                else
                    ratio = term_num / term_den;
                end
                
                PE3(t) = PE3(t) - D(i,j) * ratio * (sin(th(t,i) - th(t,j)) - sin(ths(i) - ths(j)));
            end
        end
        
        % Total Potential Energy
        PE(t) = PE1(t) + PE2(t) + PE3(t);
    end
end