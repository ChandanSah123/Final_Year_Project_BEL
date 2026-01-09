clc; clear; close all;
% 1. Setup Data
global Yint_post;
load('data1.mat'); 
load('Y_all.mat');
load('CCT_TimeDomain.mat', 'CCT_TD');
load('Fault_Info.mat');
%CCT_TD=CCT_TD+1;
T = data(:,51);
npts=length(T);
num_bus=39;
num_gen=10;
g=num_gen;
idx_ang = 1:10;       % Columns for Angle
idx_spd = 11:20;      % Columns for Speed
idx_pm  = 21:30;      % Columns for Mech Power
idx_pe  = 31:40;

% Parameters
H = [42.0;      30.3;   35.8;     28.6;   26.0;     34.8;  26.4;       24.3;   34.5;       500.0]; % Example
M = 2 * H / (2 * pi * 60); 
M_tot = sum(M);
Ws=(2*pi*60);
% Extract Data
delta = data(:, idx_ang)*(pi/180);
omega = data(:, idx_spd)*Ws; % Ensure this is speed DEVIATION (w - 1.0) or (w - w0)
Pm = data(:, idx_pm);
Pe = data(:, idx_pe);
%figure;plot(delta);
%figure; plot (Pe);

%% 2. Efficient COI Calculation (Vectorized)
% Calculate Center of Inertia (COI) without loops
d_COI = (delta * M) / M_tot; % Matrix multiplication (N x 10) * (10 x 1) -> (N x 1)
w_COI = (omega * M) / M_tot;

% Transform to COI frame
theta = delta - d_COI;
w_tilde = omega - w_COI;
w = w_tilde;
th=theta;
%figure;plot(theta);
%% Solving for post fault SEP using Load Flow
idx_gen = 30:39;
idx_load=setdiff(1:num_bus,idx_gen);
v_gen=[1.0475 0.9820 0.9831 0.9972 1.0123 1.0493 1.0635 1.0278 1.0265 1.0300]';
slack_bus=31;
xd_p=     [0.025, 0.05, 0.045, 0.035, 0.089, 0.04, 0.044, 0.045, 0.045, 0.004]';
Pgen=zeros(num_bus,1);
Pgen(idx_gen)=Pm(1,1:num_gen);
YN_post=Y_post;
x_eq_post=NR_ss(YN_post,Pgen,idx_load,v_gen,slack_bus);
V_eq_post=x_eq_post(num_bus+1:end).*(cos(x_eq_post(1:num_bus))+1i*sin(x_eq_post(1:num_bus)));
I_eq_post=YN_post*V_eq_post;
Pgen_post=real(V_eq_post(idx_gen).*conj(I_eq_post(idx_gen)));
Eeq_post=abs(V_eq_post(idx_gen)+1i*xd_p.*I_eq_post(idx_gen));
delta_eq_post=angle(V_eq_post(idx_gen)+1i*xd_p.*I_eq_post(idx_gen));
x_eq_post=[delta_eq_post-M'*delta_eq_post/M_tot; zeros(num_gen,1)];
ths=x_eq_post(1:num_gen);
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
%% 3. DIRECT CUEP & MOD IDENTIFICATION (Using CCT Snapshot)
fault_start_time = 1.0; 
t_cct_absolute1 = fault_start_time + CCT_TD+0.2; %slightly higher than cct
t_cct_absolute=fault_start_time + CCT_TD;

fprintf('Fault Duration: %.4fs | Absolute Clearing Time: %.4fs\n', CCT_TD, t_cct_absolute1);

[~, idx_cct] = min(abs(T - t_cct_absolute));

fprintf('Snapshot taken at simulation step: %d (t = %.4fs)\n', idx_cct, T(idx_cct));

theta_guess = theta(idx_cct, :)';
[sorted_angles, sort_idx] = sort(theta_guess, 'descend');
% We assume 'w_tilde' is your speed variable in COI frame
w_guess = w_tilde(idx_cct, :)'; 
sorted_speeds = w_guess(sort_idx);

% 3. Construct MOD_sort_data 
% Format: [Gen_Index, Angle_Val, Speed_Val]
MOD_sort_data = [sort_idx, sorted_angles, sorted_speeds];

fprintf('Identified MOD Candidates (Sorted by Angle):\n');
fprintf('Gen: %d | Ang: %.4f | Spd: %.4f\n', MOD_sort_data');

% We calculate the Power Injection (Pi) first
Pi = Pgen(1:num_gen) - (real(diag(Y1))) .* ((E(:)).^2);
VPE = Calculate_PE(npts, g, Pi, C, D, th, ths);
%% 4. ITERATIVE MOD SEARCH
results_table = []; 
best_error = 9999;

% Initialize "Best" variables
best_MOD_group = [];
best_CCT_TEF = 0;
best_theta_cuep = [];
best_Vcr = 0;

sorted_gen_indices = MOD_sort_data(:, 1); 

for k = 1:num_gen
    fprintf('\n--- Testing MOD Candidate: Top %d Generator(s) ---\n', k);
    
    % A. Calculate theta_u (Corner Point)
    theta_u = Calculate_theta_u(k, MOD_sort_data, num_gen, ths, H);
    current_MOD_indices = MOD_sort_data(1:k, 1); 
    
    % 0. Safety Check
    if k == num_gen
        fprintf('   Skipping k=%d (Full system cannot separate).\n', k);
        continue; 
    end

    % B. Optimization: Robust CUEP Solver (Layered Approach)
    
    % 1. Define Inputs
    Pm_vec = Pm(1, 1:num_gen)';   
    E_vec  = Eeq_post;            
    M_vec  = M;                   
    
    % 2. Function Handle (Anchored Mismatch)
    cuep_fun = @(th) cuep_mismatch_anchored(th, Pm_vec, E_vec, Y1, M_vec, M_tot);

    % 3. STRATEGY 1: lsqnonlin (Fast, Gradient-based)
    options_lsq = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
                               'Display', 'off', 'FunctionTolerance', 1e-6, ...
                               'StepTolerance', 1e-6, 'MaxFunctionEvaluations', 1000);
                           
    % Start from Corner Point (theta_u)
    fprintf('   Solving CUEP (Level 1: lsqnonlin)... ');
    [theta_sol, resnorm, ~, exitflag] = lsqnonlin(cuep_fun, theta_u, [], [], options_lsq);

    % 4. STRATEGY 2: fmincon SQP (The "Nuclear Option")
    % If Strategy 1 failed (high residual), we minimize the error square sum directly.
    if resnorm > 0.1
        fprintf('\n      > Level 1 stuck (Res: %.2f). Engaging Level 2 (fmincon SQP)...', resnorm);
        
        % Objective: Minimize sum of squared mismatches
        obj_fun = @(th) sum(cuep_mismatch_anchored(th, Pm_vec, E_vec, Y1, M_vec, M_tot).^2);
        
        % Constraints: Bound angles within +/- 60 degrees of the corner point to prevent drift
        lb = theta_u - pi/3;
        ub = theta_u + pi/3;
        
        options_sqp = optimoptions('fmincon', 'Algorithm', 'sqp', ...
                                   'Display', 'off', 'MaxFunctionEvaluations', 2000, ...
                                   'OptimalityTolerance', 1e-6);
        
        try
            [theta_sqp, fval_sqp] = fmincon(obj_fun, theta_u, [], [], [], [], lb, ub, [], options_sqp);
            
            % If SQP found a better point, take it
            if fval_sqp < resnorm
                theta_sol = theta_sqp;
                resnorm = fval_sqp;
            end
        catch
            fprintf(' (Optimization Error)');
        end
    end
    
    % 5. Final Processing & SKIPPING LOGIC
    % Re-center results to COI one last time to be precise
    d_offset = sum(M_vec .* theta_sol) / M_tot;
    theta_cuep_mod = theta_sol - d_offset;
    
    % Final Residual Check (Power mismatch only)
    final_mismatch_vec = cuep_mismatch_anchored(theta_cuep_mod, Pm_vec, E_vec, Y1, M_vec, M_tot);
    real_residual = sum(abs(final_mismatch_vec(1:num_gen))); 
    
    if real_residual < 0.5
        fprintf(' Converged. Res: %.2e\n', real_residual);
    else
        % --- THE MODIFICATION IS HERE ---
        fprintf(' \n      > Solvers Failed (Res: %.2e). Skipping this MOD.\n', real_residual);
        continue; % Jump to next 'k' immediately, skipping Vcr calculation
    end
    
    % C. Transient Energy & CCT (Only runs if solver converged)
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
   
    % D. SAVE BEST CANDIDATE
    if err < best_error
        best_error = err;
        best_MOD_group = current_MOD_indices;
        best_CCT_TEF = cct_tef;
        best_theta_cuep = theta_cuep_mod; 
        best_Vcr = Vcr_candidate;        
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
new_entry.Fault_Location =Tripped_Line; % Adjust this based on your simulation setup
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

