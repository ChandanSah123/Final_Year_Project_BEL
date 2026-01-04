clc; clear; close all;
% 1. Setup Data
global Yint_post;
load('data1.mat'); 
load('Y_all.mat');
load('CCT_TimeDomain.mat', 'CCT_TD')
%CCT_TD=CCT_TD+1;
T = data(:,16);
npts=length(T);
num_bus=9;
num_gen=3;
g=num_gen;
idx_ang = 1:3;       % Columns for Angle
idx_spd = 4:6;      % Columns for Speed
idx_pm  = 7:9;      % Columns for Mech Power
idx_pe  = 10:12;      % Columns for Elec Power

% Parameters
H = [23.64;6.4;3.01]; % Example
M = 2 * H / (2 * pi * 60); 
M_tot = sum(M);
Ws=(2*pi*60);
% Extract Data
delta = data(:, idx_ang)*(pi/180);
omega = data(:, idx_spd)*Ws; % Ensure this is speed DEVIATION (w - 1.0) or (w - w0)
Pm = data(:, idx_pm);
Pe = data(:, idx_pe);

%% 2. Efficient COI Calculation (Vectorized)
d_COI = (delta * M) / M_tot; % Matrix multiplication (N x 10) * (10 x 1) -> (N x 1)
w_COI = (omega * M) / M_tot;

% Transform to COI frame
theta = delta - d_COI;
w_tilde = omega - w_COI;
w = w_tilde;
th=theta;
%% Solving for post fault SEP using Load Flow
idx_gen=[1 2 3];
idx_load=setdiff(1:num_bus,idx_gen);
v_gen=[1.04;1.025;1.025];
slack_bus=1;
xd_p=[0.0608;0.1198;0.1813];
Pgen=zeros(9,1);
Pgen(1:3)=Pm(1,1:3);
YN_post=Y_post;
x_eq_post=NR_ss(YN_post,Pgen,idx_load,v_gen,slack_bus);
V_eq_post=x_eq_post(num_bus+1:end).*(cos(x_eq_post(1:num_bus))+1i*sin(x_eq_post(1:num_bus)));
I_eq_post=YN_post*V_eq_post;
Pgen_post=real(V_eq_post(idx_gen).*conj(I_eq_post(idx_gen)));
Eeq_post=abs(V_eq_post(idx_gen)+1i*xd_p.*I_eq_post(idx_gen));
delta_eq_post=angle(V_eq_post(idx_gen)+1i*xd_p.*I_eq_post(idx_gen));
x_eq_post=[delta_eq_post-M'*delta_eq_post/M_tot; zeros(num_gen,1)];
ths=x_eq_post(1:num_gen);
%ths=[-0.1782 0.5309 0.2711]';
%% krons reduced matrix
 Y1=Yint_post;
E=Eeq_post;
 %E=[1.054 1.050 1.017];
 C=zeros(g,g);
 D=zeros(g,g);
 for i=1:g
    for j=1:g
         C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j)); 
    end
 end
%% DIRECT CUEP & MOD IDENTIFICATION (Using CCT Snapshot)
fault_start_time = 1.0; 
t_cct_absolute1 = fault_start_time + CCT_TD+0.2; %slightly higher than cct
t_cct_absolute=fault_start_time + CCT_TD;
fprintf('Fault Duration: %.4fs | Absolute Clearing Time: %.4fs\n', CCT_TD, t_cct_absolute1);
[~, idx_cct] = min(abs(T - t_cct_absolute));
fprintf('Snapshot taken at simulation step: %d (t = %.4fs)\n', idx_cct, T(idx_cct));
theta_guess = theta(idx_cct, :)';
[sorted_angles, sort_idx] = sort(theta_guess, 'descend');
% 2. Extract the corresponding Speeds for those sorted generators
w_guess = w_tilde(idx_cct, :)'; 
sorted_speeds = w_guess(sort_idx);
MOD_sort_data = [sort_idx, sorted_angles, sorted_speeds];
fprintf('Identified MOD Candidates (Sorted by Angle):\n');
fprintf('Gen: %d | Ang: %.4f | Spd: %.4f\n', MOD_sort_data');

% We calculate the Power Injection (Pi) first
Pi = Pgen(1:num_gen) - (real(diag(Y1))) .* ((E(:)).^2);



%%
VPE = Calculate_PE(npts, g, Pi, C, D, th, ths);
%% 4. ITERATIVE MOD SEARCH
results_table = []; 
best_error = 9999;

% Initialize "Best" variables to ensure they exist even if loop fails
best_MOD_group = [];
best_CCT_TEF = 0;
best_theta_cuep = [];
best_Vcr = 0;

sorted_gen_indices = MOD_sort_data(:, 1); 

for k = 1:num_gen
    fprintf('\n--- Testing MOD Candidate: Top %d Generator(s) ---\n', k);
    mod_indx = sorted_gen_indices(1:k); 
    
    % A. Calculate theta_u
    theta_u = Calculate_theta_u1(k, MOD_sort_data, num_gen, ths, H);
    current_MOD_indices = MOD_sort_data(1:k, 1); 

    % B. Optimization (SQP)
    Pm_row = Pm(1, 1:num_gen);   
    E_row  = Eeq_post.';         
    H_col  = H(:);               

    % Center Guess & Constraints
    x0 = theta_u - (sum(theta_u .* H_col) / sum(H_col));
    Aeq = H_col'; beq = 0;            
    lb = -4*pi * ones(g, 1); ub = 4*pi * ones(g, 1);
    
    opts = optimset('Algorithm', 'sqp', 'Display', 'off', ...
                    'MaxFunEvals', 50000, 'MaxIter', 2000, ...
                    'TolCon', 1e-10, 'TolX', 1e-10, 'TolFun', 1e-10);
    
    nonlin_con = @(x) SEPfunction(x, Pm_row, E_row, C, D, H, Y1);
    obj_fun = @(x) 0;

    fprintf('   Solving for CUEP...\n');
    [x_sol, fval, exitflag, output] = fmincon(obj_fun, x0, [], [], Aeq, beq, lb, ub, nonlin_con, opts);

    if exitflag <= 0
        warning('   Optimization failed. Skipping candidate.');
        continue;
    end
    
    % Post-Process CUEP
    theta_cuep_mod = x_sol - (sum(x_sol .* H_col) / sum(H_col));
    [~, mism] = nonlin_con(theta_cuep_mod);
    fprintf('   CUEP Found. Max Mismatch: %.2e p.u.\n', max(abs(mism)));

    % C. Transient Energy & CCT
    Vcr_candidate = Calculate_PE_single_point(theta_cuep_mod, ths, Pi, C, D, g);
    [KE, KE_corr_candidate] = Calculate_KE(npts, g, H, Ws, w, current_MOD_indices);
    
    V_total = VPE + KE_corr_candidate;
    delV_candidate = Vcr_candidate - V_total;

    tcr_idx_candidate = find(delV_candidate < 0);
    
    if isempty(tcr_idx_candidate)
        fprintf('   System Stable (V < Vcr). No CCT found.\n');
        cct_tef = NaN;
        err = 9999;
    else
        tcr_idx = tcr_idx_candidate(1);
        cct_tef = T(tcr_idx);
        err = abs(cct_tef - t_cct_absolute);
        fprintf('   CCT_TEF: %.4fs | CCT_TD: %.4fs | Error: %.4f\n', cct_tef, t_cct_absolute, err);
    end
    
    results_table = [results_table; k, cct_tef, err, Vcr_candidate];
   
    % --- D. SAVE BEST CANDIDATE ---
    % We must save the parameters of the BEST iteration to the database
    if err < best_error
        best_error = err;
        best_MOD_group = current_MOD_indices;
        best_CCT_TEF = cct_tef;
        best_theta_cuep = theta_cuep_mod; % Vital: Save the correct CUEP
        best_Vcr = Vcr_candidate;         % Vital: Save the correct Vcr
    end
end

%% 7. BUILD/UPDATE OFFLINE LOOKUP TABLE (DATABASE)
% Matches structure in paper: [Fault Location | KE Signature | MOD | Vcr]

db_file = 'Offline_Database.mat';

% 1. Calculate KE Signature at Fault Clearing Time
% The signature is the Normalized Kinetic Energy of each machine at t_clear
t_clearing_time = fault_start_time + CCT_TD; 
[~, idx_clear] = min(abs(T - t_clearing_time));

w_at_clearing = w(idx_clear, :)'; % Speed deviation at clearing
KE_raw = 0.5 * M .* (w_at_clearing.^2); % Individual KE
KE_total = sum(KE_raw);

% Normalized Signature (This is what you match against in real-time)
KE_Signature = KE_raw / KE_total; 

% 2. Define the Entry Structure
new_entry.Fault_Location = 'Bus 9'; % Adjust this based on your simulation setup
new_entry.KE_Signature   = KE_Signature; % 3x1 Vector (The "Fingerprint")
new_entry.MOD_Generators = best_MOD_group;
new_entry.Critical_Energy= best_Vcr;
new_entry.CUEP_Angles    = best_theta_cuep;
new_entry.CCT_TEF        = best_CCT_TEF;

% 3. Display Data to be Saved
fprintf('\n--- Saving Best Result to Database ---\n');
fprintf('Fault Loc: %s\n', new_entry.Fault_Location);
fprintf('Best MOD: Gen [ %s]\n', num2str(best_MOD_group'));
fprintf('KE Signature: [ %.4f  %.4f  %.4f ]\n', KE_Signature');
fprintf('Critical Energy (Vcr): %.4f\n', best_Vcr);

% 4. Load & Append
if exist(db_file, 'file')
    load(db_file, 'TEF_Database');
    
    % Optional: Check for duplicates (Simple check based on Fault Loc)
    % Ideally, you would check if this specific case already exists.
    TEF_Database = [TEF_Database; new_entry];
else
    TEF_Database = new_entry;
end

save(db_file, 'TEF_Database');
fprintf('Entry added. Total Entries in DB: %d\n', length(TEF_Database));