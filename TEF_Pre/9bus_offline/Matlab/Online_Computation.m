clc; clear; close all;

%% 1. SETUP & LOADING
% -------------------------------------------------------------------------
% Load the same system data and the Offline Database we just built
global Yint_post;
load('data1.mat');              % "Real-Time" Data Stream
load('Y_all.mat');              % Admittance Matrices
load('Offline_Database.mat');   % The Look-up Table
load('Fault_Info.mat');         % To know which fault occurred (Bus/Line)

fprintf('--- ONLINE TRANSIENT STABILITY ASSESSMENT ---\n');

% System Parameters (Must match Offline Config)
num_bus = 9;
num_gen = 3;
g = num_gen;
H = [23.64; 6.4; 3.01]; 
M = 2 * H / (2 * pi * 60); 
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

%% 2. CALCULATE POST-FAULT SEP (ONLINE)
% -------------------------------------------------------------------------
% In a real setting, we run a fast load flow or use a state estimator
% to find the stable equilibrium point (SEP) for the post-fault topology.
fprintf('1. Calculating Post-Fault Equilibrium (SEP)...\n');

idx_gen = [1 2 3];
idx_load = setdiff(1:num_bus, idx_gen);
v_gen = [1.04; 1.025; 1.025];
slack_bus = 1;
xd_p = [0.0608; 0.1198; 0.1813];

Pgen = zeros(9,1);
Pgen(1:3) = Pm(1, 1:3); % Assume Pm constant
YN_post = Y_post;

% Run Newton-Raphson
x_eq_post = NR_ss(YN_post, Pgen, idx_load, v_gen, slack_bus);
V_eq_post = x_eq_post(num_bus+1:end) .* (cos(x_eq_post(1:num_bus)) + 1i*sin(x_eq_post(1:num_bus)));
I_eq_post = YN_post * V_eq_post;

% Calculate Internal SEP Angles (ths)
delta_eq_post = angle(V_eq_post(idx_gen) + 1i*xd_p.*I_eq_post(idx_gen));
ths = delta_eq_post - (M' * delta_eq_post / M_tot); % COI Reference

% Calculate Y-bus Reduced Parameters (C and D matrices)
Y1 = Yint_post;
E = abs(V_eq_post(idx_gen) + 1i*xd_p.*I_eq_post(idx_gen));
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
% -------------------------------------------------------------------------
% We will test the system state at a specific time "t_online".
% For testing, let's pick the Clearing Time from the simulation data.
% (In real deployment, this block comes from PMU stream)

% Let's test at t = Fault Clearing Time (Critical Moment)
fault_start_time = 1.0;
% Assuming we want to assess stability at the moment the fault clears:
% We look at the fault info to define the clearing time used in simulation
current_clearing_time = 1.0 + 0.2330; % Example: Test near the CCT found earlier
% Or pick a time from the T vector
[~, idx_online] = min(abs(T - current_clearing_time));

t_online = T(idx_online);
fprintf('2. Acquiring Data Snapshot at t = %.4fs\n', t_online);

% Get Measurements at this instant
delta_online = delta(idx_online, :);
omega_online = omega(idx_online, :);

% Convert to COI Frame
d_COI_online = (delta_online * M) / M_tot;
w_COI_online = (omega_online * M) / M_tot;

theta_online = (delta_online - d_COI_online)'; % Transpose to (3x1)
w_tilde_online = (omega_online - w_COI_online)'; 


%% 4. SIGNATURE GENERATION & KNN MATCHING (K > 1)
% -------------------------------------------------------------------------
fprintf('3. Generating Online KE Signature & Running KNN Search...\n');

% A. Generate Online Signature
KE_raw_online = 0.5 * M .* (w_tilde_online.^2);
KE_total_online = sum(KE_raw_online);
KE_Sig_Online = KE_raw_online / KE_total_online;

fprintf('   Online Signature: [ %.4f  %.4f  %.4f ]\n', KE_Sig_Online');

% --- KNN SETUP ---
K = 5; % Number of neighbors to check (Adjust as per paper, usually 5 or 10)
distances = [];     % Store [Distance, Database_Index]
candidates = [];    % Store valid entries

% B. Calculate Distances for ALL Entries
for i = 1:length(TEF_Database)
    entry = TEF_Database(i);
    
    % Filter by Fault Location
    if strcmp(entry.Fault_Location, Fault_Location)
        % Euclidean Distance
        dist = norm(entry.KE_Signature - KE_Sig_Online);
        
        % Store: [Distance, Index in DB]
        distances = [distances; dist, i];
    end
end

if isempty(distances)
    error('No matching fault location found in Database!');
end

% C. Sort and Select Top K
% Sort by distance (Ascending)
sorted_distances = sortrows(distances, 1);

% Select Top K (Handle case where DB has fewer than K entries)
num_available = size(sorted_distances, 1);
K_actual = min(K, num_available);
top_K_indices = sorted_distances(1:K_actual, 2);

fprintf('   Top %d Neighbors found:\n', K_actual);

% D. Voting / Selection Mechanism
% We need to find which MOD Group appears most frequently in the top K
vote_list = {}; 
vote_counts = [];

for k = 1:K_actual
    db_idx = top_K_indices(k);
    this_mod = TEF_Database(db_idx).MOD_Generators;
    
    % Convert MOD array to string for easy counting (e.g., "[2 3]")
    mod_str = num2str(this_mod(:)'); 
    
    fprintf('     Rank %d (Dist %.4f): MOD [%s]\n', k, sorted_distances(k,1), mod_str);
    
    % Check if we already have this MOD in our vote list
    found = false;
    for v = 1:length(vote_list)
        if strcmp(vote_list{v}, mod_str)
            vote_counts(v) = vote_counts(v) + 1;
            found = true;
            break;
        end
    end
    
    if ~found
        vote_list{end+1} = mod_str;
        vote_counts(end+1) = 1;
    end
end

% Find Winner
[max_votes, winner_idx] = max(vote_counts);
winning_mod_str = vote_list{winner_idx};

% Retrieve the Best Match (Closest neighbor that belongs to the winning group)
% We prefer the closest distance that *agrees* with the majority vote.
best_match_idx = -1;
for k = 1:K_actual
    db_idx = top_K_indices(k);
    this_mod = TEF_Database(db_idx).MOD_Generators;
    mod_str = num2str(this_mod(:)');
    
    if strcmp(mod_str, winning_mod_str)
        best_match_idx = db_idx;
        break; % Found the closest one that matches the majority
    end
end

Matched_Entry = TEF_Database(best_match_idx);
fprintf('   \n   WINNING MOD: [%s] with %d/%d votes.\n', winning_mod_str, max_votes, K_actual);
fprintf('   Using parameters from DB Entry #%d (Closest of the winners).\n', best_match_idx);
fprintf('   Retrieved Vcr: %.4f\n', Matched_Entry.Critical_Energy);


%% 5. STABILITY ASSESSMENT (TEF CALCULATION)
% -------------------------------------------------------------------------
fprintf('4. Calculating Transient Energy Margin...\n');

% A. Calculate Potential Energy (V_PE) at Current State
% Using the same function as offline, but with current angles theta_online
V_PE_current = Calculate_PE_single_point(theta_online, ths, Pi, C, D, g);

% B. Calculate Corrected Kinetic Energy (V_KE_corr)
% Using the MOD retrieved from the database
MOD_gens = Matched_Entry.MOD_Generators;
M_eq = sum(M(MOD_gens));
w_eq = sum(M(MOD_gens) .* w_tilde_online(MOD_gens)) / M_eq;
V_KE_corr_current = 0.5 * M_eq * (w_eq^2);

% C. Total Energy
V_total_current = V_PE_current + V_KE_corr_current;

% D. Energy Margin
Energy_Margin = Matched_Entry.Critical_Energy - V_total_current;

fprintf('\n   --- RESULTS ---\n');
fprintf('   Current V_PE      : %.4f\n', V_PE_current);
fprintf('   Current V_KE_corr : %.4f\n', V_KE_corr_current);
fprintf('   Total Energy (V)  : %.4f\n', V_total_current);
fprintf('   Critical Energy   : %.4f\n', Matched_Entry.Critical_Energy);
fprintf('   ---------------------------\n');
fprintf('   ENERGY MARGIN     : %.4f\n', Energy_Margin);

if Energy_Margin > 0
    fprintf('   PREDICTION        : STABLE \n');
else
    fprintf('   PREDICTION        : UNSTABLE \n');
end

%% Helper Function (Ensure this is in your path or at bottom of script)
function V_pe = Calculate_PE_single_point(theta, ths, Pi, C, D, g)
    % Calculates PE for a single time step (theta is g x 1 vector)
    
    term1 = 0;
    for i = 1:g
        term1 = term1 - Pi(i) * (theta(i) - ths(i));
    end
    
    term2 = 0;
    for i = 1:g-1
        for j = i+1:g
            % Energy stored in transmission lines (relative to SEP)
            % Integral of -[Cij*sin(th_ij) + Dij*cos(th_ij)] d(theta)
            
            % At State (theta)
            val_curr = -C(i,j)*cos(theta(i)-theta(j)) + D(i,j)*sin(theta(i)-theta(j));
            
            % At SEP (ths)
            val_sep  = -C(i,j)*cos(ths(i)-ths(j)) + D(i,j)*sin(ths(i)-ths(j));
            
            % Integration path implies: -(Val_curr - Val_sep)
            term2 = term2 - (val_curr - val_sep);
        end
    end
    
    V_pe = term1 + term2;
end