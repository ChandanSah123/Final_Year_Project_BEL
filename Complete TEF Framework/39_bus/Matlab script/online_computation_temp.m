clc; clear; close all;
tic
%% 1. SETUP & LOADING
load('data1.mat'); 
load('Y_all.mat');
load('CCT_TimeDomain.mat', 'CCT_TD');
load('Offline_Database.mat', 'TEF_Database');
load('Tcl.mat', 't_cl');
load('Fault_Info.mat');
 load('shedable_gen.mat');
%CCT_TD=CCT_TD+1;
num_bus=39;
num_gen=10;
g=num_gen;
H = [42.0;      30.3;   35.8;     28.6;   26.0;     34.8;  26.4;       24.3;   34.5;       500.0]; % Example
M = 2 * H / (2 * pi * 60); 
E= [1.0862; 1.1197; 1.1228; 1.0675; 1.2504; 1.2023; 1.0534; 1.0584; 1.1017; 1.0344]';
xd_p= [0.025, 0.05, 0.045, 0.035, 0.089, 0.04, 0.044, 0.045, 0.045, 0.004]';
gen_buses = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39];
load_buses = setdiff(1:num_bus, gen_buses);
load_adm_matrix = [
3,3.0593-0.0228i;
4,4.9995-1.8398i;
7,2.3636-0.8492i;
8,5.2871-1.7826i;
12,0.0754-0.885i;
15,3.1574-1.5096i;
16,3.1529-0.3067i;
18,1.5054-0.2858i;
20,6.419-1.0528i;
21,2.6431-1.1093i;
23,2.3978-0.8196i;
24,2.9336+0.8746i;
25,2.0098-0.4235i;
26,1.2648-0.1547i;
27,2.638-0.7088i;
28,1.875-0.2512i;
29,2.5782-0.2446i;
31,0.0954-0.0477i;% Load on Generator Bus
39,10.4063-2.3565i % Load on Generator Bus
];
Y1 = Yint_post;
Ybus_post=Y_post;
M_tot = sum(M);
Ws=(2*pi*60);
T = data(:,51);
npts=length(T);
idx_ang = 1:10;       % Columns for Angle
idx_spd = 11:20;      % Columns for Speed
idx_pm  = 21:30;      % Columns for Mech Power
idx_pe  = 31:40;

% Parameters


% Extract Data
delta = data(:, idx_ang)*(pi/180);
omega = data(:, idx_spd)*Ws; % Ensure this is speed DEVIATION (w - 1.0) or (w - w0)
Pm = data(:, idx_pm);
Pe = data(:, idx_pe);
%% 3. COI REFERENCE FRAME TRANSFORMATION
d_COI = (delta * M) / M_tot;
w_COI = (omega * M) / M_tot;

theta = delta - d_COI;
w_tilde = omega - w_COI;
w = w_tilde;
th = theta;

%% 4. CALCULATE POST-FAULT SEP 
fprintf('1. Calculating Post-Fault Equilibrium (SEP)...\n');
Pgen = zeros(num_gen,1);
Pgen(1:num_gen) = Pm(1, 1:num_gen); 

% Calculate Y-bus Reduced Parameters (C and D matrices)
C = zeros(g,g);
D = zeros(g,g);
for i = 1:g
    for j = 1:g
        C(i,j) = E(i)*E(j)*imag(Y1(i,j));
        D(i,j) = E(i)*E(j)*real(Y1(i,j)); 
    end
end
Pi = Pgen(1:num_gen) - (real(diag(Y1))) .* ((E(:)).^2);

% Finding SEP using fmincon
Pm_row = Pm(1, 1:num_gen);   
obj_fun_sep = @(x) 1; 
x0_sep = theta(100,:);

opts_sep = optimset('Algorithm', 'sqp', 'Display', 'off');
nonlin_con_sep = @(x) SEPfunction(x, Pm_row, E, C, D, H, Y1);

[x_sol1, ~, ~, ~] = fmincon(obj_fun_sep, x0_sep, [], [], [], [], -pi*ones(g,1), pi*ones(g,1), nonlin_con_sep, opts_sep);
ths = x_sol1 - (sum(x_sol1(:).* H(:)) / sum(H));

%% 5. REAL-TIME DATA ACQUISITION & kNN MOD MATCHING
fault_start_time = 1.0;
data_time_online = fault_start_time + 0.1; 
[~, idx_online_snap] = min(abs(T - data_time_online));
t_online = T(idx_online_snap);

fprintf('2. Acquiring Data Snapshot at t = %.4fs\n', t_online);
delta_online = theta(idx_online_snap, :);
omega_online = w(idx_online_snap, :);

KE_online = 0.5 .* M' .* (omega_online).^2;
KE_tot = sum(KE_online);
KE_Norm = KE_online / KE_tot;

% kNN Search against offline database
X = vertcat([TEF_Database.KE_Signature].');
Y_knn = KE_Norm;
[idx_knn, Euc_D] = knnsearch(X, Y_knn, "K", 5, "Distance", "euclidean");
Sorted_Mod = {TEF_Database(idx_knn).MOD_Generators};

%% 6. FINDING CUEP VECTOR FOR RANKED MODS (Bug Fixed Here)
CUEP_Candidate = zeros(length(Sorted_Mod), num_gen);
H_col = H(:);               

fprintf('Solving for CUEP with COI Constraint...\n');
for k = 1:length(Sorted_Mod)
    MOD = Sorted_Mod{k};
    theta_u = Calculate_theta_u_online(MOD, num_gen, ths, H);
    
    % --- BUG FIX: Variable 'x' is now correctly passed to Calculate_Fsum ---
    obj_fun_cuep = @(x) Calculate_Fsum(x, g, Pi, C, D, H);
    x0_cuep = theta_u;
    
    Aeq = H';       
    beq = 0; 
    lb = -2*pi * ones(g, 1); ub =  2*pi * ones(g, 1);
    
    opts_cuep = optimset('Algorithm', 'interior-point', 'Display', 'off');
    
    [x_sol_cuep, ~, exitflag_cuep, ~] = fmincon(obj_fun_cuep, x0_cuep, [], [], Aeq, beq, lb, ub, [], opts_cuep);
    theta_cuep_mod = x_sol_cuep - (sum(x_sol_cuep .* H_col) / sum(H_col));
    
    if exitflag_cuep <= 0
        warning('CUEP Optimization did not converge. ExitFlag: %d', exitflag_cuep);
    else
        % Calculate actual mismatch to verify finding a true saddle point
        mismatch_val = Calculate_Fsum(theta_cuep_mod, g, Pi, C, D, H);
        fprintf('CUEP Found. Power Mismatch (Sum of Squares): %e\n', mismatch_val);
    end
    
    CUEP_Candidate(k,:) = theta_cuep_mod;
end

%% 7. MOD IDENTIFICATION & MARGIN ASSESSMENT
faultclearing_time = 1 + t_cl;
[~, idx_clearing] = min(abs(T - faultclearing_time));
wcl = w(idx_clearing, :);
thetacl = theta(idx_clearing, :);

limit = min(length(Sorted_Mod), 3);
Vcr_candidate = zeros(limit, 1);
KE_corr_candidate = zeros(limit, 1);
V_total = zeros(limit, 1);
del_VPE_candidate = zeros(limit, 1);
norm_delVPE = zeros(limit, 1);

VPE = Calculate_PE_single_point(thetacl, ths, Pi, C, D, g);

for i = 1:limit
    MOD = Sorted_Mod{i};
    Vcr_candidate(i) = Calculate_PE_single_point(CUEP_Candidate(i,:), ths, Pi, C, D, g);
    KE_corr_candidate(i) = Calculate_KE_online(g, H, Ws, wcl, MOD);
    V_total(i) = VPE + KE_corr_candidate(i);
    
    del_VPE_candidate(i) = Vcr_candidate(i) - V_total(i);
    norm_delVPE(i) = del_VPE_candidate(i) / KE_corr_candidate(i);
end

ResultTable = table((1:limit)', Vcr_candidate, KE_corr_candidate, del_VPE_candidate, norm_delVPE, ...
    'VariableNames', {'Rank', 'V_CUEP', 'KE_Corr', 'Margin_DeltaV', 'Norm_Margin'});
disp(ResultTable);

[~, idx_mod] = min(norm_delVPE);
Correct_MOD = Sorted_Mod{idx_mod};

% Assessment of First Swing Stability
thcuep = CUEP_Candidate(idx_mod,:);
Vcr = Calculate_PE_single_point(thcuep, ths, Pi, C, D, g);
VPEcl = Calculate_PE_single_point(thetacl, ths, Pi, C, D, g);
KE_corr_cl = Calculate_KE_online(g, H, Ws, wcl, Correct_MOD);
Vcl = VPEcl + KE_corr_cl;
del_V = Vcr - Vcl;

fprintf("Energy margin = %f\n", del_V);
if del_V < 0
    Shed = 1;
    fprintf("Energy margin is negative. System is Unstable. Taking control action...\n");
else
    Shed = 0;
    fprintf("Energy margin is positive. System is Stable. No action required.\n");
end

%% 8. AUTOMATED GENERATOR SHEDDING ALGORITHM
if Shed == 1
    
    threshold = 0.001; 
    Psh_max = 1;
    max_iter = 100;
    
    pred_data_time = faultclearing_time+0.2;
    [~, idx_online] = min(abs(T - pred_data_time));
    delta_abs = delta(idx_online, :); 
    omega_abs = omega(idx_online, :);
    display(delta_abs);
    display(omega_abs);
    
    % --- INITIALIZE WARM START VARIABLES ---
    ths_guess = ths(:);
    thcuep_guess = CUEP_Candidate(idx_mod,:)';
    
    k = 1;
    Psh(k) = 0;
    VPE_val(k) = Calculate_PE_single_point(theta(idx_online,:), ths, Pi, C, D, g);
    VKE_val(k) = Calculate_KE_online(g, H, Ws, w(idx_online,:), Correct_MOD);
    Vtot(k) = VPE_val(k) + VKE_val(k);
    Vsh(k) = Vcr - Vtot(k);
    
    fprintf('Iter %d: Shed %05.2f%% -> Vcr: %.4f, Vtot: %.4f, Margin: %.4f\n', ...
            k, Psh(k)*100, Vcr, Vtot(k), Vsh(k));


 % DETERMINATION OF SHEDDABLE GENERATORS
    data_time = 1 + t_cl;
    [~, idx_online] = min(abs(T - data_time));
    wpred = omega(idx_online,:);
    deltapred = delta(idx_online,:);
    Pgen = Pm(100, 1:num_gen);
    
    % Calculate ratio of Kinetic Energy to Mechanical Power
    KE = 0.5 * M' .* (wpred).^2;
    Ratio = KE ./ Pgen;
    [max_val, Gen_sh] = max(Ratio);  %Identify the most severely accelerated generator
    Gen_idx = Gen_sh;
   
    k = k + 1;
    Psh(k) = 0.01; % Initial 1% step
    while k < max_iter
        % A. Update Network Parameters
        H1 = H; Pgen1 = Pgen; xd_p1 = xd_p;
        
        H1(Gen_idx) = H(Gen_idx) * (1 - Psh(k));
        Pgen1(Gen_idx) = Pgen(Gen_idx) * (1 - Psh(k));
        xd_p1(Gen_idx) = xd_p(Gen_idx) / (1 - Psh(k));
        
        Y1_shed = calculate_kron_reduction(Ybus_post, xd_p1, gen_buses, load_adm_matrix);
        
        C_new = zeros(g,g); D_new = zeros(g,g);
        for i=1:g
            for j=1:g
                C_new(i,j) = E(i)*E(j)*imag(Y1_shed(i,j));
                D_new(i,j) = E(i)*E(j)*real(Y1_shed(i,j)); 
            end
        end
        Pi1 = Pgen1(:) - (real(diag(Y1_shed))) .* ((E(:)).^2);
        
        % B. UPDATE SEP (Using Warm Start)
        obj_fun_sep = @(x) 1;
        nonlin_con_sep = @(x) SEPfunction(x, Pgen1(1:num_gen)', E, C_new, D_new, H1, Y1_shed);
        opts_opt = optimset('Algorithm','sqp','Display','off');
        
        [x_sep, ~, exit1] = fmincon(obj_fun_sep, ths_guess, [], [], [], [], -pi*ones(g,1), pi*ones(g,1), nonlin_con_sep, opts_opt);
        if exit1 > 0
            ths_new = x_sep - (sum(x_sep(:).* H1(:)) / sum(H1));
            ths_guess = ths_new; % WARM START: Save for next iter
        else
            ths_new = ths_guess; 
        end

        % C. UPDATE CUEP (Using Warm Start and Fixed Objective Function)
        obj_fun_cuep_shed = @(x) Calculate_Fsum(x, g, Pi1, C_new, D_new, H1);
        opts_opt = optimset('Algorithm','interior-point','Display','off');
        [x_cuep, ~, exit2] = fmincon(obj_fun_cuep_shed, thcuep_guess, [], [], H1', 0, -2*pi*ones(g,1), 2*pi*ones(g,1), [], opts_opt);
        
        if exit2 > 0
             thcuep_new = x_cuep - (sum(x_cuep(:).* H1(:)) / sum(H1));
             Vcr_new = Calculate_PE_single_point(thcuep_new, ths_new, Pi1, C_new, D_new, g);
             thcuep_guess = thcuep_new; % WARM START: Save for next iter
        else
             Vcr_new = Vcr; 
        end
        
        % D. RE-CALCULATE COI FRAME (Mass Changed)
        M1 = 2 * H1 / Ws;
        M_tot1 = sum(M1);
        
        d_COI_new = (delta_abs * M1) / M_tot1;
        w_COI_new = (omega_abs * M1) / M_tot1;
        
        deltapred_new = delta_abs - d_COI_new;
        wpred_new     = omega_abs - w_COI_new;
        
        % E. CALCULATE MARGIN
        VPE_val(k) = Calculate_PE_single_point(deltapred_new, ths_new, Pi1, C_new, D_new, g);
        VKE_val(k) = Calculate_KE_online(g, H1, Ws, wpred_new, Correct_MOD);
        Vtot(k) = VPE_val(k) + VKE_val(k);
        Vsh(k) = Vcr_new - Vtot(k); 
        
        fprintf('Iter %d: Shed %05.2f%% -> Vcr: %.4f, Vtot: %.4f, Margin: %.4f\n', ...
            k, Psh(k)*100, Vcr_new, Vtot(k), Vsh(k));
            
 % F. CHECK STABILITY AND CALCULATE NEXT SHEDDING AMOUNT
        
        if Vsh(k) >= threshold  
            % FLOWCHART DIAMOND 1: "Vsh[k] < threshold ?" -> YES
            fprintf('>>> Analytical stable margin found at %.2f%% shedding.\n', Psh(k)*100);
            
            % FLOWCHART DIAMOND 2: "Psh[k] < Psh_MAX ?"
            if Psh(k) < Psh_max
                % YES PATH: "Shed Psh amount of generation" -> "Command to Unit"
                fprintf('>>> Executing Physical Shed Command: %.2f%% of Gen %d.\n', Psh(k)*100, Gen_idx);
                break; % End loop, command sent successfully
            else
                % NO PATH: "Take other control actions"
                warning('>>> Stable margin found, but required shed (%.2f%%) exceeds Psh_MAX limit. Taking other control actions.', Psh(k)*100);
                break; % End loop, trigger secondary safety systems
            end
            
        else 
             if k > max_iter || Psh(k) > 2.0 
                 warning('Secant method diverged. Aborting.');
                 break;
             end
              % (Optional safeguard to prevent MATLAB from freezing if secant diverges forever, 
            % acting as an emergency loop-breaker, though technically not in the flowchart)

            % FLOWCHART DIAMOND 1: "Vsh[k] < threshold ?" -> NO
            
            % FLOWCHART BLOCK: "Interpolate/Extrapolate Psh value at V=0 i.e. P1 using two points nearest to V=0"
            if k > 1
                [~, sorted_indices] = sort(abs(Vsh(1:k)));
                idx1 = sorted_indices(1);
                idx2 = sorted_indices(2);
                
                if abs(Psh(idx1) - Psh(idx2)) < 1e-4
                    P_new = Psh(k) + 0.05; 
                else
                    slope = (Vsh(idx1) - Vsh(idx2)) / (Psh(idx1) - Psh(idx2));
                    if abs(slope) < 1e-4
                        slope = sign(Vsh(k)) * 0.01; 
                    end
                    target_V = 0; % Root finding target
                    P_new = Psh(idx1) + (target_V - Vsh(idx1)) / slope;
                end
                
                % Numeric bounds to prevent wild secant extrapolations
                delta_P = P_new - Psh(k);
                delta_P = max(min(delta_P, 0.10), -0.10); 
                P_new = Psh(k) + delta_P;
            else
                P_new = Psh(k) + 0.02; % Step size if only 1 point exists
            end
            
            % Prevent physically impossible negative shedding guesses
            P_new = max(P_new, 0); 
            
            % FLOWCHART BLOCK: "k = k + 1"
            k = k + 1;
            
            % Set next Psh value for the beginning of the next loop iteration 
            % (Aligns with Flowchart: "Psh[k] = P1")
            Psh(k) = P_new;
            
        end
    end
end

toc