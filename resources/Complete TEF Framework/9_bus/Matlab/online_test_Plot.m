clc; clear; close all;
tic
%% 1. SETUP & LOADING
load('data1.mat');              % "Real-Time" Data Stream
load('Y_all.mat');     
load('CCT_TimeDomain.mat', 'CCT_TD');  % Admittance Matrices
load('Offline_Database.mat', 'TEF_Database')  % The Look-up Table
load('Fault_Info.mat');   % To know which fault occurred (Bus/Line)
load('Tcl.mat');
%load('shedable_gen.mat');

fprintf('--- ONLINE TRANSIENT STABILITY ASSESSMENT ---\n');
num_bus = 9;
num_gen = 3;
g = num_gen;
H = [23.64; 6.4; 3.01]; 
M = 2 * H / (2 * pi * 60); 
E=[1.0587 1.0503 1.0170];
xd_p = [0.0608; 0.1198; 0.1813];
gen_buses = [1 2 3];
load_buses = setdiff(1:num_bus, gen_buses);
load_adm_matrix = [
    5, (1.26 - 0.504i);
    6, (0.877 - 0.292i);
    8, (0.969 - 0.339i)
];
Y1 = Yint_post;
Ybus_post=Y_post;
M_tot = sum(M);
Ws = (2 * pi * 60);
% Extract Data
T = data(:,16);
idx_ang = 1:3;
idx_spd = 4:6;
idx_pm  = 7:9;
idx_pe  = 10:12;

delta = data(:, idx_ang) * (pi/180);
omega = data(:, idx_spd) * Ws; 
Pm    = data(:, idx_pm);
%% Determination of the shedable generators
data_time=1+t_cl;
[~, idx_online] = min(abs(T - data_time));
wpred=omega(idx_online,:);
deltapred=delta(idx_online,:);
Pgen=Pm(100,1:num_gen);
KE=0.5*M'.*(wpred).^2;
Ratio= KE ./ Pgen;
[max_val, idx] = max(Ratio);
Gen_sh=idx;
%% COI reference
d_COI = (delta * M) / M_tot;
w_COI = (omega * M) / M_tot;
theta = delta - d_COI;
w_tilde = omega - w_COI;
w = w_tilde;
th=theta;

%% 2. CALCULATE POST-FAULT SEP 
% to find the stable equilibrium point (SEP) for the post-fault topology.
fprintf('1. Calculating Post-Fault Equilibrium (SEP)...\n');

Pgen = zeros(num_gen,1);
Pgen(1:num_gen) = Pm(1, 1:num_gen); % Assume Pm constant

%%  Calculate Y-bus Reduced Parameters (C and D matrices)
C = zeros(g,g);
D = zeros(g,g);
for i = 1:g
    for j = 1:g
        C(i,j) = E(i)*E(j)*imag(Y1(i,j));
        D(i,j) = E(i)*E(j)*real(Y1(i,j)); 
    end
end
% Calculate Power Injection (Pi)
Pi = Pgen(1:num_gen) - (real(diag(Y1))) .* ((E(:)).^2);

%% 3. SIMULATE "REAL-TIME" DATA ACQUISITION
fault_start_time = 1.0;
data_time = fault_start_time+ 0.1; % Example: Test near the CCT found earlier

[~, idx_online] = min(abs(T - data_time));

t_online = T(idx_online);
fprintf('2. Acquiring Data Snapshot at t = %.4fs\n', t_online);

% Get Measurements at this instant
delta_online = theta(idx_online, :);
omega_online = w(idx_online, :);

KE=0.5.*M'.*(omega_online).^2;
KE_tot=sum(KE);
KE_Norm=KE/KE_tot;

%% kNN search
X =vertcat([TEF_Database.KE_Signature].');
Y = KE_Norm;
[idx_knn,Euc_D]=knnsearch(X,Y,"K",5,"Distance","euclidean");

Fault = {TEF_Database(idx_knn).Fault_Location}.';
Sorted_Mod={TEF_Database(idx_knn).MOD_Generators};
%% Alternative method to find the stable equilibrium point
Pm = Pm(1, 1:num_gen);   % 1x3 Row Vector
obj_fun = @(x) 1;
x0 = theta(100,:);
%linear Constraint
A = []; 
b = []; 
Aeq = [];       
beq = []; 
%bounds
lb = -pi * ones(g, 1); ub =  pi * ones(g, 1);
%opts = optimset('Algorithm', 'interior-point');
opts = optimset('Algorithm', 'sqp','display','off');
% 5. Constraint Function Handle
nonlin_con =@(x) SEPfunction(x, Pm, E, C, D, H, Y1);
% 6. Run Optimization
fprintf('Solving for CUEP with COI Constraint...\n');
[x_sol1, fval1, exitflag1, output1] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, nonlin_con, opts);
ths = x_sol1-(sum(x_sol1(:).* H(:)) / sum(H));
%% Finding CUEP vector For each ranked MOD
CUEP_Candidate=zeros(length(Sorted_Mod),num_gen);
for k=1:length(Sorted_Mod)
MOD=Sorted_Mod{k};
theta_u = Calculate_theta_u_online(MOD, num_gen, ths, H);
Pm_row = Pm(1, 1:num_gen);   % 1x3 Row Vector
E_row  = E;         % 1x3 Row Vector
H_col  = H(:);               % 3x1 Column Vector
% objective function
obj_fun = @(x)Calculate_Fsum(x, g, Pi, C, D, H);
%initial condition
x0 = theta_u;
%linear Constraint
A = []; 
b = []; 
Aeq = H';       
beq = 0; 

%bounds

lb = -2*pi * ones(g, 1);  ub =  2*pi * ones(g, 1);


%opts = optimset('Algorithm', 'interior-point');
opts = optimset('Algorithm', 'sqp');
% 5. Constraint Function Handle
nonlin_con =[];

% 6. Run Optimization
fprintf('Solving for CUEP with COI Constraint...\n');
[x_sol, fval, exitflag, output] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, nonlin_con, opts);
theta_cuep_mod= x_sol- (sum(x_sol .* H_col) / sum(H_col));

%%
if exitflag <= 0
    warning('CUEP Optimization did not converge. ExitFlag: %d', exitflag);
else
    %[~, mism] = nonlin_con(theta_cuep_mod);
    fprintf('CUEP Found. Max Power Mismatch\n');
end
  

fprintf('Saving the controlling Unstable Equillibrium Point \n');
CUEP_Candidate(k,:)=theta_cuep_mod;
end
%% Correct MOD Identification
%Checking for few MODs
 faultclearing_time=1+t_cl;
 [~, idx_clearing] = min(abs(T - faultclearing_time));
 wcl=w(idx_clearing, :);
 thetacl=theta(idx_clearing, :);
num_candidates = length(Sorted_Mod);
limit = min(num_candidates, 3);
Vcr_candidate = zeros(limit, 1);
KE_corr_candidate = zeros(limit, 1);
V_total = zeros(limit, 1);
del_VPE_candidate = zeros(limit, 1);
norm_delVPE = zeros(limit, 1);
VPE= Calculate_PE_single_point(thetacl, ths, Pi, C, D, g);
for i=1:limit
MOD=Sorted_Mod{i};
Vcr_candidate(i) = Calculate_PE_single_point(CUEP_Candidate(i,:), ths, Pi, C, D, g);
KE_corr_candidate(i) = Calculate_KE_online(g, H, Ws, wcl, MOD);
V_total(i) = VPE + KE_corr_candidate(i);
del_VPE_candidate(i)=Vcr_candidate(i)-V_total(i);
norm_delVPE(i)=del_VPE_candidate(i)/KE_corr_candidate(i);
end

ResultTable = table((1:limit)', Vcr_candidate, KE_corr_candidate, del_VPE_candidate, norm_delVPE, ...
    'VariableNames', {'Rank', 'V_CUEP', 'KE_Corr', 'Margin_DeltaV', 'Norm_Margin'});

disp(ResultTable);
% Find the lowers normalized potential energy
[min_maargin, idx_mod]=min(norm_delVPE);
Correct_MOD=Sorted_Mod{idx_mod};
% [2 3]
%[2 3]
%% Assessment of the First Swing Stability
 thcuep=CUEP_Candidate(idx_mod,:);
 Vcr= Calculate_PE_single_point(thcuep, ths, Pi, C, D, g);
 VPEcl= Calculate_PE_single_point(thetacl, ths, Pi, C, D, g);
 KE_corr_cl = Calculate_KE_online(g, H, Ws, wcl,Correct_MOD);
 Vcl=VPEcl+KE_corr_cl;
 del_V=Vcr-Vcl;
 fprintf("Energy margin=%d\n",del_V);
if del_V<0
    Shed=1;
    fprintf("Energy margin is negative System is Unstable take control action \n");
else
    Shed=0;
    fprintf("Energy margin is positive System is Stable no Need to take control Action\n");
   
end
%% Shedding Calculation
if Shed == 1
    Gen_idx = Gen_sh;
    threshold = 0; 
    Psh_max = 1;
    max_iter = 30;
    pred_data_time = faultclearing_time+0.05;
    [~, idx_online] = min(abs(T - pred_data_time));
    delta_abs = delta(idx_online, :); 
    omega_abs = omega(idx_online, :);
    wpred=w(idx_online,:);
    thetapred=theta(idx_online,:);
    % Store the Original CUEP for the solver's initial guess
 
    k = 1;
    Psh(k) = 0;
    % Initial Margin Calculation
    VPE(k) = Calculate_PE_single_point(thetapred, ths, Pi, C, D, g);
    VKE(k) = Calculate_KE_online(g, H, Ws, wpred, Correct_MOD);
    Vtot(k) = VPE(k) + VKE(k);
    Vsh(k) = Vcr - Vtot(k);
     fprintf('Iter %d: Shed %.2f%% -> Vcr: %.2f, Vtot: %.2f, Margin: %.4f\n', ...
            k, Psh(k)*100, Vcr, Vtot(k), Vsh(k));
    % Start Iteration
    k = k + 1;
    Pmin = 0.01; 
    Psh(k) = Pmin;
    
    optimal_shedding_found = false;
    
   while k < max_iter
        % --- A. UPDATE NETWORK PARAMETERS ---
        H1 = H; Pgen1 = Pgen; xd_p1 = xd_p;
        
        % Apply Shedding
        H1(Gen_idx) = H(Gen_idx) * (1 - Psh(k));
        Pgen1(Gen_idx) = Pgen(Gen_idx) * (1 - Psh(k));
        xd_p1(Gen_idx) = xd_p(Gen_idx) / (1 - Psh(k));
        
        % Update Y-Bus
        Y1_shed = calculate_kron_reduction(Ybus_post, xd_p1, gen_buses, load_adm_matrix);
        
        E_new = E; 
        C_new = zeros(g,g); D_new = zeros(g,g);
        for i=1:g
            for j=1:g
                C_new(i,j) = E_new(i)*E_new(j)*imag(Y1_shed(i,j));
                D_new(i,j) = E_new(i)*E_new(j)*real(Y1_shed(i,j)); 
            end
        end
        Pi1 = Pgen1(1:num_gen) - (real(diag(Y1_shed))) .* ((E_new(:)).^2);
        
        % --- B. RE-CALCULATE COI FRAME (Using new Inertia) ---
        M1 = 2 * H1 / Ws;
        M_tot1 = sum(M1);
        
        d_COI_new = (delta_abs * M1) / M_tot1;
        w_COI_new = (omega_abs * M1) / M_tot1;
        
        deltapred_new = delta_abs - d_COI_new;
        wpred_new     = omega_abs - w_COI_new;

        % --- C. CALCULATE ENERGY AFTER SHEDDING (As per flowchart) ---
        % Notice we use the ORIGINAL 'ths' and 'Vcr' here, just updated grid matrices
        VPE(k) = Calculate_PE_single_point(deltapred_new, ths, Pi1, C_new, D_new, g);
        VKE(k) = Calculate_KE_online(g, H1, Ws, wpred_new, Correct_MOD);
        Vtot(k) = VPE(k) + VKE(k);
        
        % Calculate Margin using the ORIGINAL Vcr calculated before the loop
        Vsh(k) = Vcr - Vtot(k); 
        
        fprintf('Iter %d: Shed %.2f%% -> Vcr: %.2f, Vtot: %.2f, Margin: %.4f\n', ...
            k, Psh(k)*100, Vcr, Vtot(k), Vsh(k));

        % --- D. DECISION & INTERPOLATION (Matching the Flowchart) ---
        if Vsh(k) >= threshold  % Flowchart: Vsh[k] < threshold (adjusted for your positive margin logic)
            final_cmd = Psh(k);
            fprintf('>>> Stable margin found at %.2f%% shedding.\n', Psh(k)*100);
            fprintf('>>> Command to Unit: Shed %.2f%% of gen %d.\n', Psh(k)*100, Gen_sh);
            break;
        else 
            if Psh(k) >= Psh_max
                 warning('Required shedding exceeds limit. Take Other Action.');
                 break;
            end
            
            % Flowchart: "Interpolate/Extrapolate Psh value at V=0 using two points"
            if k > 1
                y2 = Vsh(k);     x2 = Psh(k);
                y1 = Vsh(k-1);   x1 = Psh(k-1);
                
                slope = (y2 - y1) / (x2 - x1);
                
                if abs(slope) < 1e-4
                    slope = sign(y2) * 0.01; 
                end
                
                delta_P = (threshold - y2) / slope;
                
                P_new = x2 + delta_P;
            else
                % If k=2 (first shedding pass), we don't have two shed points to interpolate yet
                P_new = Psh(k) + 0.02; 
            end
            
            P_new = max(min(P_new, Psh_max), 0);
            if P_new <= Psh(k) && Vsh(k) < 0
                 P_new = Psh(k) + 0.01; 
            end
            
            k = k + 1;
            Psh(k) = P_new;
        end
   end
end
