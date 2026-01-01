function F = PowerMismatch(th, Pi, C, D, nGen)
    F = zeros(nGen, 1);
    for i = 1:nGen
        Pe_i = 0;
        for j = 1:nGen
            % Summing C*sin(d) + D*cos(d)
            Pe_i = Pe_i + C(i,j)*sin(th(i)-th(j)) + D(i,j)*cos(th(i)-th(j));
        end
        F(i) = Pi(i) - Pe_i; 
    end
end