clc; clear; close all;
%% Setup 
load('data1.mat'); 
load('Y_all.mat');
load('CCT_TimeDomain.mat', 'CCT_TD');
load('Fault_Info.mat');
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
E=[1.0587 1.0503 1.0170];
xd_p=[0.0608;0.1198;0.1813];
M_tot = sum(M);
Ws=(2*pi*60);
% Extract Data
delta = data(:, idx_ang)*(pi/180);
omega = data(:, idx_spd)*Ws; % Ensure this is speed DEVIATION (w - 1.0) or (w - w0)
Pm = data(:, idx_pm);
Pe = data(:, idx_pe);

%% 2. Efficient COI Calculation (Vectorized)
d_COI = (delta * M) / M_tot;
w_COI = (omega * M) / M_tot;
% Transform to COI frame
theta = delta - d_COI;
w_tilde = omega - w_COI;
w = w_tilde;
th=theta;
%% Solving for post fault SEP using Load Flow
Pgen=zeros(num_gen,1);
Pgen(1:num_gen)=Pm(1,1:num_gen);
%% krons reduced matrix
Y1=Yint_post;
 C=zeros(g,g);
 D=zeros(g,g);
 for i=1:g
    for j=1:g
         C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j)); 
    end
 end
 Pi = Pgen(1:num_gen) - (real(diag(Y1))) .* ((E(:)).^2);
 %% Alternative method to find the stable equilibrium point
Pm = Pm(1, 1:num_gen);   % 1x3 Row Vector
E=[1.0587 1.0503 1.0170];
H = H(:); 
obj_fun = @(x) 1;
x0 = theta(100,:);
%linear Constraint
A = []; 
b = []; 
Aeq = [];       
beq = []; 
lb = -pi * ones(g, 1); 
ub =  pi * ones(g, 1);
opts = optimset('Algorithm', 'sqp');
%opts = optimset('Algorithm', 'interior-point');
% 5. Constraint Function Handle
nonlin_con =@(x) SEPfunction(x, Pm, E, C, D, H, Y1);
% 6. Run Optimization
fprintf('Solving for CUEP with COI Constraint...\n');
[x_sol1, fval1, exitflag1, output1] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, nonlin_con, opts);
ths = x_sol1-(sum(x_sol1(:) .* H(:)) / sum(H));

%% DIRECT CUEP & MOD IDENTIFICATION (Using CCT Snapshot)
fault_start_time = 1.0; 
t_cct_absolute1 = fault_start_time + CCT_TD+0.1; %slightly higher than cct
t_cct_absolute=fault_start_time + CCT_TD;
fprintf('Fault Duration: %.4fs | Absolute Clearing Time: %.4fs\n', CCT_TD, t_cct_absolute1);
[~, idx_cct] = min(abs(T - t_cct_absolute1));
fprintf('Snapshot taken at simulation step: %d (t = %.4fs)\n', idx_cct, T(idx_cct));
theta_guess = theta(idx_cct, :)';
[sorted_angles, sort_idx] = sort(theta_guess, 'descend');
% 2. Extract the corresponding Speeds for those sorted generators
w_guess = w_tilde(idx_cct, :)'; 
sorted_speeds = w_guess(sort_idx);
MOD_sort_data = [sort_idx, sorted_angles, sorted_speeds];
fprintf('Identified MOD Candidates (Sorted by Angle):\n');
fprintf('Gen: %d | Ang: %.4f | Spd: %.4f\n', MOD_sort_data');

%%
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

for k = 1:num_gen-1
    fprintf('\n--- Testing MOD Candidate: Top %d Generator(s) ---\n', k);
    mod_indx = sorted_gen_indices(1:k); 
    
    % A. Calculate theta_u (Corner Point)
    theta_u = Calculate_theta_u(k, MOD_sort_data, num_gen, ths, H);
    current_MOD_indices = MOD_sort_data(1:k, 1); 
    
% 1. Setup Inputs
% Ensure inputs match the dimensions required by SEPfunction (Row Vectors)
Pm_row = Pm(1, 1:num_gen);   % 1x3 Row Vector
E_row  = E;                  % 1x3 Row Vector
H_col  = H(:);               % 3x1 Column Vector
% objective function
obj_fun = @(x)Calculate_Fsum(x, g, Pi, C, D, H);
%initial condition
x0 = theta_u;
%linear Constraint
A = []; 
b = []; 
Aeq = [];       
beq = []; 

%bounds

lb = -2*pi * ones(g, 1); 
ub =  2*pi * ones(g, 1);


opts = optimset('Algorithm', 'sqp');
%opts = optimset('Diplay','iter','Algorithm','sqp');
% 5. Constraint Function Handle
nonlin_con =[];

% 6. Run Optimization
fprintf('Solving for CUEP with COI Constraint...\n');
[x_sol, fval, exitflag, output] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, nonlin_con, opts);

% 7. Post-Process
theta_cuep_mod = x_sol - (sum(x_sol(:) .* H_col(:)) / sum(H_col));

% 8. Convergence Check
if exitflag <= 0
    warning('CUEP Optimization did not converge. ExitFlag: %d', exitflag);
else
   % [~, mism] = nonlin_con(theta_cuep_mod);
   % fprintf('CUEP Found. Max Power Mismatch: %.2e p.u.\n', max(abs(mism)));
end

 % C. Transient Energy & CCT
 % C. Transient Energy & CCT
    Vcr_candidate = Calculate_PE_single_point(theta_cuep_mod, ths, Pi, C, D, g);
    [KE, KE_corr_candidate] = Calculate_KE(npts, g, H, Ws, w, current_MOD_indices);
    
    % Calculate Energy Margin Trajectory
    V_total = VPE + KE_corr_candidate;
    delV_candidate = Vcr_candidate - V_total; % Stable if delV > 0
    
    % Find first index where margin becomes negative
    unstable_indices = find(delV_candidate < 0);
   if isempty(unstable_indices)
    fprintf('System is stable for the entire duration of T.\n');
    cct_tef=NaN;
    err=9999;
    
else
    idx_unstable = unstable_indices(1);
    if idx_unstable <= 1
         cct_tef = 0;
         fprintf('System is unstable from the start.');
       
    else
        t_stable   = T(idx_unstable - 1);
        t_unstable = T(idx_unstable);
        margin_stable   = delV_candidate(idx_unstable - 1);
        margin_unstable = delV_candidate(idx_unstable);
        % Linear Interpolation formula to find where Margin = 0
        slope = (margin_unstable - margin_stable) / (t_unstable - t_stable);
        dt = (0 - margin_stable) / slope;
        cct_tef = t_stable + dt;
    end
     err = abs(cct_tef - t_cct_absolute);
     fprintf('CCT Found: %.4fs (Interpolated)\n', cct_tef);
     fprintf('   Vcr: %.4f | KE_corr (at fault): %.4f | VPE: %.4f \n', Vcr_candidate, KE_corr_candidate(idx_cct),VPE(idx_cct));
     fprintf('   CCT_TEF: %.5fs | CCT_TD: %.4fs | Error: %.4f\n', cct_tef, t_cct_absolute, err);
  end
    
    results_table = [results_table; k, cct_tef, err, Vcr_candidate];
   
    % --- D. SAVE BEST CANDIDATE ---
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
% The signature is the Normalized Kinetic Energy of each machine at 0.1
% second of fault
t_clearing_time = fault_start_time + 0.1; 
[~, idx_clear] = min(abs(T - t_clearing_time));

w_at_clearing = w(idx_clear, :)'; % Speed deviation at clearing
KE_raw = 0.5 * M .* (w_at_clearing.^2); % Individual KE
KE_total = sum(KE_raw);

% Normalized Signature (This is what you match against in real-time)
KE_Signature = KE_raw / KE_total; 
% 3. Define Entry Structure
new_entry.Fault_Location    = Tripped_Line; 
new_entry.MOD_Generators    = best_MOD_group; 
new_entry.KE_Signature      = KE_Signature;
new_entry.Vcr               = best_Vcr;
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


 
 