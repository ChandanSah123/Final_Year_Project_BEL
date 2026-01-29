function [c, ceq] = SEPfunction(x, Pm, E, C, D, H, Y)

 
    c = [];
    g = numel(x); % Automatically detect system size (e.g., 3 or 10)
    Pcoi = 0;
    for i = 1:g 
        for j = i+1:g
            Pcoi = Pcoi + D(i,j)*cos(x(i)-x(j));
        end
    end
    Pcoi = Pcoi * 2;
    
    ceq = zeros(g, 1);
    
    for i = 1:g
        % Now this loop is safe because g matches the size of x
        for k = 1:g
            if x(k) > pi
                x(k) = -(2*pi - x(k));
            end
            if x(k) < -pi
                x(k) = (2*pi + x(k));
            end
        end
        
        for j = 1:g        
            if i ~= j
                ceq(i) = ceq(i) - (C(i,j)*sin(x(i)-x(j)) + D(i,j)*cos(x(i)-x(j)));
            else
                ceq(i) = ceq(i) + Pm(i) - E(i)^2 * real(Y(i,i));
            end
        end
        
        % Ensure explicit row/column handling for the COI term
        % We use 'sum(Pm)' assuming Pm is a row vector from your Main script fixes
        term_loss = sum(Pm - E.^2 .* real(diag(Y))'); 
        ceq(i) = ceq(i) - (H(i)/sum(H)) * term_loss + (H(i)/sum(H)) * Pcoi;
    end
end