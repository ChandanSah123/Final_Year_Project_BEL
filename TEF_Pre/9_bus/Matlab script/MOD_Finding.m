clc; clear; close all;
% 1. Setup Data
global Yint_post;
load('data1.mat'); 
load('Y_all.mat');
load('CCT_TimeDomain.mat', 'CCT_TD')
CCT_TD=CCT_TD+1;
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
% Calculate Center of Inertia (COI) without loops
d_COI = (delta * M) / M_tot; % Matrix multiplication (N x 10) * (10 x 1) -> (N x 1)
w_COI = (omega * M) / M_tot;

% Transform to COI frame
theta = delta - d_COI;
w_tilde = omega - w_COI;
w = w_tilde;
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
%% Fth formulation
ths=x_eq_post(1:3);
th=theta;
Pi=Pgen(1:3)-(real(diag(Y1))).*((E(:)).^2);
 fth = Calculate_fth(ths, th, Pi, npts, g, C, D, H);
% figure;plot(fth,'g');
 plot(fth);

%% Detect PEBS Crossing (Modified)
% 1. DEFINE START POINT
fault_start_time = 1.0;  % Change this to match your exact fault inception time
t_sim=T;
[MOD_sort_data, t_pebs_idx, PEBS_Crossing] = detect_pebs_crossing(fth, t_sim, fault_start_time, npts, th, w_tilde);

%% performing the gradinet decent/sigma F(th) minimization using pebs data
tic
thpebs = PEBS_Crossing(:,2); % 10x1 vector
options = odeset('MaxStep',0.01,'InitialStep',0.01);

% PASS PARAMETERS: @(t,y) Integrateth(t, y, Pi, C, D, H)
[Tp, Y] = ode45(@(t,y) Integrateth(t, y, Pi, C, D, H), [0 9], thpebs, options);

yp = Y;
toc
theta_cuep_pebs = yp(end, :)';  
fprintf('CUEP Angles (radians) BCU: \n');
disp(theta_cuep_pebs);


%% 5. MOD Identification Algorithm (Iterative Lookup)
fprintf('\n==============================================\n');
fprintf('     STARTING MOD IDENTIFICATION LOOP\n');
fprintf('==============================================\n');

% Pre-calculate VPE (Potential Energy) as it does not depend on MOD
VPE = Calculate_PE(npts, g, Pi, C, D, th, ths);

% Prepare to store results
results_table = []; % Columns: [Num_Machines, CCT_TEF, Error, Vcr]
best_error = 9999;
best_MOD_group = [];
best_CCT_TEF = 0;

sorted_gen_indices = MOD_sort_data(:, 1); 
% LOOP: Try Top 1 machine, then Top 2 .....
for k = 1:num_gen
    fprintf('\n--- Testing MOD Candidate: Top %d Generator(s) ---\n', k);
    current_MOD_indices = sorted_gen_indices(1:k); 
    % 2. Calculate theta_u
    theta_u = Calculate_theta_u(k, MOD_sort_data, num_gen, ths, H);
    % 3. Calculate CUEP for this specific MOD (Time Domain BCU)
    options = odeset('MaxStep',0.01,'InitialStep',0.01);
    [Tm, Ym] = ode45(@(t,y) Integrateth(t, y, Pi, C, D, H), [0 20], theta_u, options);
    theta_cuep_mod = Ym(end, :)';
    % 4. Calculate Critical Energy (Vcr) for this CUEP
    Vcr_candidate = Calculate_PE_single_point(theta_cuep_mod, ths, Pi, C, D, g);
    % 5. Calculate Corrected Kinetic Energy (Depends on MOD!)
    [~, KE_corr_candidate] = Calculate_KE(npts, g, H, Ws, w, current_MOD_indices);
    % 6. Total Energy & Margin
    V_total = VPE + KE_corr_candidate;
    delV_candidate = Vcr_candidate - V_total;
    % 7. Find CCT_TEF (Zero Crossing of Margin)
    cct_tef = NaN; % Default if no crossing found
    % Start search after fault inception
    idx_search = find(T >= 1.0, 1); 
    if isempty(idx_search), idx_search = 2; end
    for t = idx_search+1 : npts
        if delV_candidate(t) < 0 && delV_candidate(t-1) >= 0
            % Linear Interpolation
            y_prev = delV_candidate(t-1);
            y_curr = delV_candidate(t);
            fraction = y_prev / (y_prev - y_curr);
            cct_tef = T(t-1) + fraction * (T(t) - T(t-1));
            break; % Found the first crossing
        end
    end
    
    % 8. Error Calculation
    if isnan(cct_tef)
        fprintf('   -> System remained Stable (V < Vcr). No CCT found.\n');
        err = 9999;
    else
        err = abs(cct_tef - CCT_TD);
        fprintf('   -> CCT_TEF: %.4fs | CCT_TD: %.4fs | Error: %.4f\n', cct_tef, CCT_TD, err);
    end
    
    % Store in Table
    results_table = [results_table; k, cct_tef, err, Vcr_candidate];
    
    % Check if this is the best one
    if err < best_error
        best_error = err;
        best_MOD_group = current_MOD_indices;
        best_CCT_TEF = cct_tef;
    end
end

%% 7. BUILD/UPDATE OFFLINE DATABASE
% This block saves the "Learned" Information into a database file.

db_file = 'Offline_Database.mat';

% 1. Calculate the Kinetic Energy "Fingerprint" (Normalized KE) for this case
% We use the state at Fault Clearing Time (approx) or Max KE point
[~, idx_clear] = min(abs(T - (1.0 + 0.22))); % Example: Fault clears at 1.22s
w_online = w(idx_clear, :)';
KE_raw = 0.5 * M .* (w_online.^2);
KE_fingerprint = KE_raw / sum(KE_raw); % Normalized (Sum = 1)

% 2. Create the Entry Structure
new_entry.Fault_Location = 'Bus 9'; % You can automate this name
new_entry.MOD_Generators = best_MOD_group;
new_entry.CUEP_Angles    = theta_cuep_mod; % The CUEP we found using that MOD
new_entry.Critical_Energy= Vcr_candidate;  % The Vcr for that MOD
new_entry.KE_Signature   = KE_fingerprint; % The Search Key

% 3. Append to Database
if exist(db_file, 'file')
    load(db_file, 'TEF_Database');
    TEF_Database = [TEF_Database; new_entry];
else
    TEF_Database = new_entry;
end

save(db_file, 'TEF_Database');
fprintf('Entry added to Offline_Database.mat. Total Entries: %d\n', length(TEF_Database));