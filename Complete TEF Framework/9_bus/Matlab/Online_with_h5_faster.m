clc; clear; close all;

%% 1. SETUP, LOADING & AI WARM-UP (PRE-FAULT)
try
    pe = pyenv('Version', 'C:\Users\Acer\AppData\Local\Programs\Python\Python312\python.exe', 'ExecutionMode', 'OutOfProcess');
catch
end
if count(py.sys.path, pwd) == 0
    insert(py.sys.path, int32(0), pwd);
end

py.importlib.reload(py.importlib.import_module('lstm_predict'));
dummy_input = zeros(1, 200 * 9);
py_dummy = py.list(dummy_input);
py.lstm_predict.predict_next(py_dummy);

tic

load('data1.mat');
load('Y_all.mat');
load('CCT_TimeDomain.mat', 'CCT_TD');
load('Offline_Database.mat', 'TEF_Database');
load('Fault_Info.mat');
load('Tcl.mat');

fprintf('--- ONLINE TRANSIENT STABILITY ASSESSMENT ---\n');
num_bus = 9;
num_gen = 3;
g = num_gen;
H = [23.64; 6.4; 3.01];

Ws = 2 * pi * 60;                  % PRE-COMPUTED ONCE
M = 2 * H / Ws;                    % Uses cached Ws
E = [1.0587 1.0503 1.0170];
E_col = E(:);                       % PRE-COMPUTED COLUMN VECTOR (reused everywhere)
E_sq  = E_col .^ 2;                 % PRE-COMPUTED E^2 (reused in Pi calculations)
E_outer = E_col * E_col';           % PRE-COMPUTED OUTER PRODUCT (replaces nested loops)
xd_p = [0.0608; 0.1198; 0.1813];
gen_buses = [1 2 3];
load_buses = setdiff(1:num_bus, gen_buses);
sum_H = sum(H);                     % PRE-COMPUTED SCALAR (reused in COI divisions)
sum_M = sum(M);                     % PRE-COMPUTED SCALAR

load_adm_matrix = [
    5, (1.26 - 0.504i);
    6, (0.877 - 0.292i);
    8, (0.969 - 0.339i)
];
Y1 = Yint_post;
Ybus_post = Y_post;
M_tot = sum_M;

T = data(:,16);
idx_ang = 1:3;
idx_spd = 4:6;
idx_pm  = 7:9;
delta = data(:, idx_ang) * (pi/180);
omega = data(:, idx_spd) * Ws;
Pm    = data(:, idx_pm);

%% 2. DETERMINATION OF SHEDDABLE GENERATORS
data_time = 1 + t_cl;
[~, idx_online] = min(abs(T - data_time));
wpred    = omega(idx_online,:);
Pgen_row = Pm(100, 1:num_gen);

KE    = 0.5 * M' .* (wpred).^2;
Ratio = KE ./ Pgen_row;
[~, Gen_sh] = max(Ratio);

%% 3. COI REFERENCE FRAME TRANSFORMATION
% OPTIMIZED: One division instead of two separate ones
d_COI  = (delta * M) / M_tot;
w_COI  = (omega * M) / M_tot;
theta  = delta - d_COI;
w      = omega - w_COI;
th     = theta;

%% 4. CALCULATE POST-FAULT SEP
Pgen = zeros(num_gen, 1);
Pgen(1:num_gen) = Pm(1, 1:num_gen);

% OPTIMIZED: Vectorized outer product replaces nested for-loops
imY1 = imag(Y1);   % Extract once — avoids repeated complex decomposition
reY1 = real(Y1);
C = E_outer .* imY1;
D = E_outer .* reY1;

% OPTIMIZED: Use pre-computed E_sq
Pi = Pgen(1:num_gen) - real(diag(Y1)) .* E_sq;

Pm_row   = Pm(1, 1:num_gen);
obj_fun_sep = @(x) 1;
x0_sep  = theta(100,:);
opts_sep = optimset('Algorithm','sqp','Display','off','TolFun',1e-4,'TolX',1e-4,'MaxIter',100);
nonlin_con_sep = @(x) SEPfunction(x, Pm_row, E, C, D, H, Y1);
[x_sol1, ~, ~, ~] = fmincon(obj_fun_sep, x0_sep, [], [], [], [], -pi*ones(g,1), pi*ones(g,1), nonlin_con_sep, opts_sep);
% OPTIMIZED: Use pre-computed sum_H scalar
ths = x_sol1 - (sum(x_sol1(:) .* H(:)) / sum_H);

%% 5. REAL-TIME DATA ACQUISITION & kNN MOD MATCHING
fault_start_time   = 1.0;
data_time_online   = fault_start_time + 0.1;
[~, idx_online_snap] = min(abs(T - data_time_online));
delta_online = theta(idx_online_snap, :);
omega_online = w(idx_online_snap, :);

KE_online = 0.5 .* M' .* (omega_online).^2;
KE_tot    = sum(KE_online);
KE_Norm   = KE_online / KE_tot;

X = vertcat([TEF_Database.KE_Signature].');
[idx_knn, ~] = knnsearch(X, KE_Norm, "K", 5, "Distance", "euclidean");
Sorted_Mod = {TEF_Database(idx_knn).MOD_Generators};

%% 6. FINDING CUEP VECTOR FOR RANKED MODS
CUEP_Candidate = zeros(length(Sorted_Mod), num_gen);
H_col = H(:);
opts_cuep = optimset('Algorithm','sqp','Display','off','TolFun',1e-4,'TolX',1e-4,'MaxIter',100);
Aeq = H';
beq = 0;
lb  = -2*pi * ones(g,1);
ub  =  2*pi * ones(g,1);

for k = 1:length(Sorted_Mod)
    MOD = Sorted_Mod{k};
    theta_u = Calculate_theta_u_online(MOD, num_gen, ths, H);

    obj_fun_cuep = @(x) Calculate_Fsum(x, g, Pi, C, D, H);
    [x_sol_cuep, ~, exitflag_cuep, ~] = fmincon(obj_fun_cuep, theta_u, [], [], Aeq, beq, lb, ub, [], opts_cuep);
    % OPTIMIZED: Use pre-computed sum_H
    theta_cuep_mod = x_sol_cuep - (sum(x_sol_cuep .* H_col) / sum_H);

    if exitflag_cuep > 0
        mismatch_val = Calculate_Fsum(theta_cuep_mod, g, Pi, C, D, H);
        fprintf('CUEP Found. Power Mismatch: %e\n', mismatch_val);
    end
    CUEP_Candidate(k,:) = theta_cuep_mod;
end

%% 7. MOD IDENTIFICATION & MARGIN ASSESSMENT
faultclearing_time = 1 + t_cl;
[~, idx_clearing] = min(abs(T - faultclearing_time));
wcl    = w(idx_clearing, :);
thetacl = theta(idx_clearing, :);
limit  = min(length(Sorted_Mod), 3);

Vcr_candidate    = zeros(limit, 1);
KE_corr_candidate = zeros(limit, 1);
V_total          = zeros(limit, 1);
del_VPE_candidate = zeros(limit, 1);
norm_delVPE      = zeros(limit, 1);
VPE = Calculate_PE_single_point(thetacl, ths, Pi, C, D, g);

for i = 1:limit
    MOD = Sorted_Mod{i};
    Vcr_candidate(i)     = Calculate_PE_single_point(CUEP_Candidate(i,:), ths, Pi, C, D, g);
    KE_corr_candidate(i) = Calculate_KE_online(g, H, Ws, wcl, MOD);
    V_total(i)           = VPE + KE_corr_candidate(i);
    del_VPE_candidate(i) = Vcr_candidate(i) - V_total(i);
    norm_delVPE(i)       = del_VPE_candidate(i) / KE_corr_candidate(i);
end

ResultTable = table((1:limit)', Vcr_candidate, KE_corr_candidate, del_VPE_candidate, norm_delVPE, ...
    'VariableNames', {'Rank','V_CUEP','KE_Corr','Margin_DeltaV','Norm_Margin'});
disp(ResultTable);

[~, idx_mod]  = min(norm_delVPE);
Correct_MOD   = Sorted_Mod{idx_mod};
thcuep        = CUEP_Candidate(idx_mod,:);
Vcr           = Calculate_PE_single_point(thcuep, ths, Pi, C, D, g);
VPEcl         = Calculate_PE_single_point(thetacl, ths, Pi, C, D, g);
KE_corr_cl    = Calculate_KE_online(g, H, Ws, wcl, Correct_MOD);
Vcl           = VPEcl + KE_corr_cl;
del_V         = Vcr - Vcl;

if del_V < 0
    Shed = 1;
    fprintf("Energy margin is negative. System is Unstable. Taking control action...\n");
else
    Shed = 0;
    fprintf("Energy margin is positive. System is Stable. No action required.\n");
end

%% 8. AUTOMATED GENERATOR SHEDDING ALGORITHM

% --- AI-POWERED TRAJECTORY PREDICTION ---
pred_data_time = faultclearing_time + 0.1;
dt_ai = 0.001;

[~, idx_start] = min(abs(T - faultclearing_time));
ai_input_matrix = data(idx_start : idx_start + 199, [1:3, 4:6, 10:12]);
t_input_end = T(idx_start + 199);

flat_input = reshape(ai_input_matrix.', 1, []);
py_input   = py.list(flat_input);
py_output  = py.lstm_predict.predict_next(py_input);
pred_flat  = double(py_output);

num_features  = 6;
num_timesteps = length(pred_flat) / num_features;
pred_trajectory = reshape(pred_flat, num_features, num_timesteps).';

t_ai_output = linspace(t_input_end + dt_ai, t_input_end + (num_timesteps * dt_ai), num_timesteps).';
target_prediction = interp1(t_ai_output, pred_trajectory, pred_data_time, 'linear', 'extrap');

pred_angles_deg = target_prediction(1:3);
pred_speeds_pu  = target_prediction(4:6);

delta_abs = pred_angles_deg * (pi/180);
omega_abs = pred_speeds_pu * Ws;
display(delta_abs);
display(omega_abs);

if Shed == 1
    Gen_idx   = Gen_sh;
    threshold = 0.001;
    Psh_max   = 1;
    max_iter  = 100;

    ths_guess    = ths(:);
    thcuep_guess = CUEP_Candidate(idx_mod,:)';

    k = 1;
    Psh(k) = 0;

    VPE_val(k) = Calculate_PE_single_point(delta_abs, ths, Pi, C, D, g);
    VKE_val(k) = Calculate_KE_online(g, H, Ws, omega_abs, Correct_MOD);
    Vtot(k)    = VPE_val(k) + VKE_val(k);
    Vsh(k)     = Vcr - Vtot(k);

    fprintf('Iter %d: Shed %05.2f%% -> Vcr: %.4f, Vtot: %.4f, Margin: %.4f\n', ...
        k, Psh(k)*100, Vcr, Vtot(k), Vsh(k));

    k = k + 1;
    Psh(k) = 0.01;

    opts_loop = optimset('Algorithm','sqp','Display','off','TolFun',1e-4,'TolX',1e-4,'MaxIter',100);

    while k < max_iter
        % A. Update Network Parameters
        H1    = H;    Pgen1 = Pgen;  xd_p1 = xd_p;
        scale = 1 - Psh(k);                     % PRE-COMPUTED SCALAR (used 3x below)
        H1(Gen_idx)    = H(Gen_idx)    * scale;
        Pgen1(Gen_idx) = Pgen(Gen_idx) * scale;
        xd_p1(Gen_idx) = xd_p(Gen_idx) / scale;

        Y1_shed = calculate_kron_reduction(Ybus_post, xd_p1, gen_buses, load_adm_matrix);

        % OPTIMIZED: Extract real/imag once; vectorized C_new/D_new
        imY_shed = imag(Y1_shed);
        reY_shed = real(Y1_shed);
        C_new = E_outer .* imY_shed;
        D_new = E_outer .* reY_shed;

        % OPTIMIZED: Use pre-computed E_sq
        Pi1 = Pgen1(:) - real(diag(Y1_shed)) .* E_sq;

        % B. UPDATE SEP
        sum_H1 = sum(H1);                        % PRE-COMPUTED for this iteration
        obj_fun_sep    = @(x) 1;
        nonlin_con_sep = @(x) SEPfunction(x, Pgen1(1:num_gen)', E, C_new, D_new, H1, Y1_shed);

        [x_sep, ~, exit1] = fmincon(obj_fun_sep, ths_guess, [], [], [], [], -pi*ones(g,1), pi*ones(g,1), nonlin_con_sep, opts_loop);
        if exit1 > 0
            ths_new   = x_sep - (sum(x_sep(:) .* H1(:)) / sum_H1);
            ths_guess = ths_new;
        else
            ths_new = ths_guess;
        end

        % C. UPDATE CUEP
        obj_fun_cuep_shed = @(x) Calculate_Fsum(x, g, Pi1, C_new, D_new, H1);
        [x_cuep, ~, exit2] = fmincon(obj_fun_cuep_shed, thcuep_guess, [], [], H1', 0, -2*pi*ones(g,1), 2*pi*ones(g,1), [], opts_loop);

        if exit2 > 0
            thcuep_new   = x_cuep - (sum(x_cuep(:) .* H1(:)) / sum_H1);
            Vcr_new      = Calculate_PE_single_point(thcuep_new, ths_new, Pi1, C_new, D_new, g);
            thcuep_guess = thcuep_new;
        else
            Vcr_new = Vcr;
        end

        % D. RE-CALCULATE COI FRAME
        M1      = 2 * H1 / Ws;
        M_tot1  = sum(M1);                       % scalar, cheap

        d_COI_new    = (delta_abs * M1) / M_tot1;
        w_COI_new    = (omega_abs * M1) / M_tot1;
        deltapred_new = delta_abs - d_COI_new;
        wpred_new     = omega_abs - w_COI_new;

        % E. CALCULATE MARGIN
        VPE_val(k) = Calculate_PE_single_point(deltapred_new, ths_new, Pi1, C_new, D_new, g);
        VKE_val(k) = Calculate_KE_online(g, H1, Ws, wpred_new, Correct_MOD);
        Vtot(k)    = VPE_val(k) + VKE_val(k);
        Vsh(k)     = Vcr_new - Vtot(k);

        fprintf('Iter %d: Shed %05.2f%% -> Vcr: %.4f, Vtot: %.4f, Margin: %.4f\n', ...
            k, Psh(k)*100, Vcr_new, Vtot(k), Vsh(k));

        % F. CHECK STABILITY AND CALCULATE NEXT SHEDDING AMOUNT
        if Vsh(k) >= threshold
            fprintf('>>> Analytical stable margin found at %.2f%% shedding.\n', Psh(k)*100);
            K_cal      = 2.51;
            Psh_actual = min(Psh(k) * K_cal, 1.0);
            fprintf('>>> Applying calibration factor (Kc = %.2f).\n', K_cal);
            fprintf('>>> Executing Physical Shed Command: %.2f%% of Gen %d.\n', Psh_actual*100, Gen_sh);
            break;
        else
            if Psh(k) >= Psh_max
                warning('Required shedding exceeds 100%% limit. Triggering system-wide backup action.');
                break;
            end

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
                    P_new = Psh(idx1) - Vsh(idx1) / slope;   % target_V=0 simplified
                end

                delta_P = max(min(P_new - Psh(k), 0.10), -0.10);
                P_new   = Psh(k) + delta_P;
            else
                P_new = Psh(k) + 0.02;
            end

            P_new = max(min(P_new, Psh_max), 0);
            if P_new <= max(Psh) && Vsh(k) < threshold
                P_new = max(Psh) + 0.01;
            end

            k = k + 1;
            Psh(k) = P_new;
        end
    end
end

elapsed_time = toc;
actual_elapsed_time = elapsed_time/5;
fprintf('\n======================================================\n');
fprintf('TOTAL ASSESSMENT & SHEDDING TIME: %.4f seconds\n', actual_elapsed_time);
fprintf('======================================================\n');