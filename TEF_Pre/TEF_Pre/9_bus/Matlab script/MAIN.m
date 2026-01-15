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
% Calculate Center of Inertia (COI) without loops
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
%% 3. DIRECT CUEP & MOD IDENTIFICATION (Using CCT Snapshot)
% Replaces "Fth formulation", "Detect PEBS", and "Gradient Descent"

% 1. DEFINE TIMING
% CCT_TD is the DURATION (e.g., 0.22s). 
% We need the ABSOLUTE TIME (e.g., 1.0 + 0.22 = 1.22s).
fault_start_time = 1.0; 
t_cct_absolute1 = fault_start_time + CCT_TD+0.2; %slightly higher than cct
t_cct_absolute=fault_start_time + CCT_TD;

fprintf('Fault Duration: %.4fs | Absolute Clearing Time: %.4fs\n', CCT_TD, t_cct_absolute1);

[~, idx_cct] = min(abs(T - t_cct_absolute));

fprintf('Snapshot taken at simulation step: %d (t = %.4fs)\n', idx_cct, T(idx_cct));

theta_guess = theta(idx_cct, :)';
[sorted_angles, sort_idx] = sort(theta_guess, 'descend');

% 2. Extract the corresponding Speeds for those sorted generators
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



%%
VPE = Calculate_PE(npts, g, Pi, C, D, th, ths);
mod_indx=2;
current_MOD_indices = MOD_sort_data(1:mod_indx,1); 

theta_u = Calculate_theta_u(mod_indx, MOD_sort_data, num_gen, ths, H);

%% code for optimizaton of theta_u goes here to obtain theta_cuep_mod

% 1. Setup Inputs
% Ensure inputs match the dimensions required by SEPfunction (Row Vectors)
Pm_row = Pm(1, 1:num_gen);   % 1x3 Row Vector
E_row  = Eeq_post.';         % 1x3 Row Vector
H_col  = H(:);               % 3x1 Column Vector

% 2. Center Initial Guess
% Ensure theta_u starts strictly at COI = 0
x0 = theta_u;

% 3. Configure Optimization Constraints
% Linear Equality Constraint: Sum(H * theta) = 0
Aeq = H_col';       
beq = 0;            

% 4. Optimization Options
obj_fun = @(x) 0; 
lb = -4*pi * ones(g, 1); 
ub =  4*pi * ones(g, 1);
A = []; b = []; 

opts = optimset('Algorithm', 'interior-point', ...
                'Display', 'iter', ...
                'MaxFunEvals', 50000, ...
                'MaxIter', 2000, ...
                'TolCon', 1e-9, ...
                'TolX', 1e-9, ...
                'TolFun', 1e-9);

% 5. Constraint Function Handle
nonlin_con = @(x) SEPfunction(x, Pm_row, E_row, C, D, H, Y1);

% 6. Run Optimization
fprintf('Solving for CUEP with COI Constraint...\n');
[x_sol, fval, exitflag, output] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, nonlin_con, opts);

% 7. Post-Process
theta_cuep_mod = x_sol - (sum(x_sol .* H_col) / sum(H_col));

% 8. Convergence Check
if exitflag <= 0
    warning('CUEP Optimization did not converge. ExitFlag: %d', exitflag);
else
    [~, mism] = nonlin_con(theta_cuep_mod);
    fprintf('CUEP Found. Max Power Mismatch: %.2e p.u.\n', max(abs(mism)));
end
% ... [Rest of your code] ...
%theta_cuep_mod = [-0.7451 2.3909 0.6487]';

Vcr_candidate = Calculate_PE_single_point(theta_cuep_mod, ths, Pi, C, D, g);
%Vcr_candidate = Calculate_PE(1, g, Pi, C, D, theta_u, ths);
[KE, KE_corr_candidate] = Calculate_KE(npts, g, H, Ws, w, current_MOD_indices);
V_total = VPE + KE_corr_candidate;
delV_candidate = Vcr_candidate - V_total;
tcr_idx_candidate=find(delV_candidate<0);
tcr_idx=tcr_idx_candidate(1);
tcr=T(tcr_idx);
display(theta_u);
display(theta_cuep_mod);
display(Vcr_candidate);
display(tcr);
display(Vcr_candidate);
figure; plot(delV_candidate);