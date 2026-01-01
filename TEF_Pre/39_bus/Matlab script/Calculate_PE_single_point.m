function Vcr = Calculate_PE_single_point(th_vec, ths, Pi, C, D, g)
    PE1 = sum(Pi .* (th_vec - ths));
    PE2 = 0; PE3 = 0;
    for i = 1:g-1
        for j = i+1:g
            th_ij = th_vec(i) - th_vec(j);
            ths_ij = ths(i) - ths(j);
            PE2 = PE2 + C(i,j) * (cos(th_ij) - cos(ths_ij));
            
            term_num = (th_vec(i) - ths(i)) + (th_vec(j) - ths(j));
            term_den = th_ij - ths_ij;
            if abs(term_den) < 1e-8, ratio = 0; else, ratio = term_num / term_den; end
            
            PE3 = PE3 - D(i,j) * ratio * (sin(th_ij) - sin(ths_ij));
        end
    end
    Vcr = PE1 + PE2 + PE3;
end