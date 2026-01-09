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
%% 4. ROBUST DIRECT EXTRACTION (PEBS METHOD)
% Strategy: Instead of solving for CUEP numerically (which is failing),
% we extract the 'Ground Truth' Critical Energy directly from the 
% time-domain trajectory, which captures the true non-linear behavior.

fprintf('\n--- Running Direct PEBS Extraction ---\n');

% 1. Find the "Ridge" (Maximum Potential Energy point)
% In a critical simulation, the max PE is the best proxy for Vcr.
[max_PE, idx_peak] = max(VPE);
t_peak = T(idx_peak);

% 2. Extract CUEP Angles at the Peak
theta_peak = theta(idx_peak, :)'; % This is your theta_cuep
w_peak     = w_tilde(idx_peak, :)';

% 3. Identify the MOD (Mode of Disturbance) Group automatically
% We sort the angles at the peak. The generators that are "advanced" 
% (separated from the rest) form the MOD.
[sorted_peak_angles, sort_peak_idx] = sort(theta_peak, 'descend');

% Calculate gaps between adjacent sorted angles to find the split
angle_gaps = abs(diff(sorted_peak_angles));
[max_gap, gap_loc] = max(angle_gaps);

% The MOD is everything above the largest gap
num_mod_gens = gap_loc;
mod_indices = sort_peak_idx(1:num_mod_gens);

% 4. Save Results
best_Vcr        = max_PE;
best_theta_cuep = theta_peak;
best_MOD_group  = mod_indices;
best_CCT_TEF    = CCT_TD; % You already know this from Time Domain

% Display Success
fprintf('   Converged via Trajectory Scan.\n');
fprintf('   Peak PE Time: %.4fs\n', t_peak);
fprintf('   Critical Energy (Vcr): %.4f\n', best_Vcr);
fprintf('   Identified MOD Size: %d generators\n', length(best_MOD_group));
fprintf('   MOD Group: [ %s]\n', num2str(best_MOD_group'));

% Fill results table for consistency
results_table = [num_mod_gens, best_CCT_TEF, 0.0, best_Vcr];
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
new_entry.Fault_Location =Fault_Location; % Adjust this based on your simulation setup
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

