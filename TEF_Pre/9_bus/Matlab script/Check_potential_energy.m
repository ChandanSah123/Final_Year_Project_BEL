
% 1. Define the two CUEP candidates
cuep_bcu = [ 9.1275    9.1395    9.2032    9.2543    9.4110    9.3376    9.5220  9.3265   35.2580 -106.1737];
cuep_mod = [10.0811   10.2110   10.2619   10.2599   10.4032   10.3539   10.5393  10.3936   34.0040 -110.5436];

% 2. Retrieve System Constants (Ensure these are available in workspace)
%    You need: ths (Post-fault SEP), Pi, C, D, num_gen
%    If 'ths' is not defined, use the starting point of your integration (x_eq_post)

% 3. Calculate Potential Energy for BCU Candidate
%    (Using the logic from the scalar function)
V_bcu = calculate_PE_single_point(cuep_bcu, ths, Pi, C, D, num_gen);

% 4. Calculate Potential Energy for MOD Candidate
V_mod = calculate_PE_single_point(cuep_mod, ths, Pi, C, D, num_gen);

% 5. Compare
fprintf('------------------------------------------------\n');
fprintf('Potential Energy at BCU CUEP: %.4f p.u.\n', V_bcu);
fprintf('Potential Energy at MOD CUEP: %.4f p.u.\n', V_mod);
fprintf('------------------------------------------------\n');

if V_bcu < V_mod
    disp('DECISION: Use BCU Result (Lower Energy Barrier).');
else
    disp('DECISION: Check angles physically. (Lower Energy is usually correct).');
end


% --- Helper Local Function (Paste at bottom of script) ---
function PE = calculate_PE_single_point(th_vec, ths, Pi, C, D, g)
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
    PE = PE1 + PE2 + PE3;
end