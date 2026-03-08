clc; clear; close all;
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
E= [1.0862; 1.1197; 1.1228; 1.0675; 1.2504; 1.2023; 1.0534; 1.0584; 1.1017; 1.0344]';
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
idx_gen=[1 2 3 4 5 6 7 8 9 10];
idx_load=setdiff(1:num_bus,idx_gen);
slack_bus=1;
xd_p=[0.0608;0.1198;0.1813];
Pgen=zeros(num_bus,1);
Pgen(1:num_gen)=Pm(1,1:num_gen);
YN_post=Y_post;
%% krons reduced matrix
Y1=Yint_post;
 %E=[1.054 1.050 1.017];
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
%E = Eeq_post.';         % 1x3 Row Vector

obj_fun = @(x) 1;
x0 = theta(1000,:);
%linear Constraint
A = []; 
b = []; 
Aeq = [];       
beq = []; 

%bounds

lb = -pi * ones(g, 1); 
ub =  pi * ones(g, 1);


opts = optimset('Algorithm', 'interior-point');
%opts = optimset('Diplay','iter','Algorithm','sqp');
% 5. Constraint Function Handle
nonlin_con =@(x) SEPfunction(x, Pm, E, C, D, H, Y1);

% 6. Run Optimization
fprintf('Solving for CUEP with COI Constraint...\n');
[x_sol1, fval1, exitflag1, output1] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, nonlin_con, opts);
ths = x_sol1-(sum(x_sol1(:).* H(:)) / sum(H));
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

%%
VPE = Calculate_PE(npts, g, Pi, C, D, th, ths);
%% 4. ITERATIVE MOD SEARCH
tic
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
E_row  = E;         % 1x3 Row Vector
H_col  = H(:);               % 3x1 Column Vector
% objective function
obj_fun = @(x)Calculate_Fsum(th, g, Pi, C, D, H);
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
toc
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
    tcr_idx_candidate = find(delV_candidate < 0);
    
    if isempty(tcr_idx_candidate)
        fprintf('   System Stable (V < Vcr). No CCT found.\n');
        cct_tef = NaN;
        err = 9999;
    else
        % --- INTERPOLATION BLOCK START ---
        idx = tcr_idx_candidate(1);
        
        if idx > 1
            % Linear Interpolation for Zero Crossing
            t1 = T(idx-1);
            t2 = T(idx);
            v1 = delV_candidate(idx-1); % Positive (Stable)
            v2 = delV_candidate(idx);   % Negative (Unstable)
            
            % Fraction of time step where crossing occurs
            frac = v1 / (v1 - v2); 
            cct_tef = t1 + frac * (t2 - t1);
        else
            cct_tef = T(idx); % Fallback if failure happens at step 1
        end
        % --- INTERPOLATION BLOCK END ---

        err = abs(cct_tef - t_cct_absolute);
        
        % Print detailed diagnostics to see the difference
        fprintf('   Vcr: %.4f | KE_corr (at fault): %.4f\n', Vcr_candidate, KE_corr_candidate(idx_cct));
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
toc
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

%%  Calculation of kinetic energy to power ratio
 
 
 