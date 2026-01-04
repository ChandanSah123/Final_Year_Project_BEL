clc; clear; close all;
% 1. Setup Data
global Yint_post;
load('data1.mat'); 
load('Y_all.mat');
load('CCT_TimeDomain.mat', 'CCT_TD')
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
%%
VPE = Calculate_PE(npts, g, Pi, C, D, th, ths);
mod_indx=1;
current_MOD_indices = MOD_sort_data(1:mod_indx,1); 

theta_u = Calculate_theta_u(mod_indx, MOD_sort_data, num_gen, ths, H);
%% code for optimizaton of theta_u goes here to obtain theta_cuep_mod

% 1. Setup Inputs (Sanitization)
% Ensure inputs match the dimensions required by SEPfunction (Row Vectors)
Pm_row = Pm(1, 1:num_gen);   % 1x10 Row Vector
E_row  = Eeq_post.';         % 1x10 Row Vector (Assuming Eeq_post is column)
H_col  = H(:);               % Ensure Inertia is a Column Vector for Aeq calculation

% 2. Center Initial Guess
% Ensure theta_u starts strictly at COI = 0
x0 = theta_u;

% 3. Configure Optimization Constraints
% CRITICAL FIX: Add Linear Equality Constraint for COI (Sum(H*theta) = 0)
% This removes the singularity/rotational invariance that was confusing the solver.
Aeq = H_col';       % Linear constraint: [H1 H2 ... Hn] * [th1; th2; ...]
beq = 0;            % Result must equal 0

% 4. Optimization Options
% Objective is zero (Feasibility Problem)
obj_fun = @(x) 0; 

% Bounds: Loose bounds to allow finding the saddle point
lb = -2*pi * ones(g, 1); 
ub =  2*pi * ones(g, 1);
A = []; b = []; 

% Solver Settings: 'sqp' is often better for singular/constrained problems,
% but 'interior-point' with Aeq set is very robust.
opts = optimset('Algorithm', 'interior-point', ...
                'Display', 'iter', ...        % Keep 'iter' to monitor convergence
                'MaxFunEvals', 50000, ...
                'MaxIter', 2000, ...
                'TolCon', 1e-9, ...           % Tighter tolerance for precision
                'TolX', 1e-9, ...
                'TolFun', 1e-9);

% 5. Constraint Function Handle
nonlin_con = @(x) SEPfunction(x, Pm_row, E_row, C, D, H, Y1);

% 6. Run Optimization
fprintf('Solving for CUEP with COI Constraint...\n');
[x_sol, fval, exitflag, output] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, nonlin_con, opts);

% 7. Post-Process
% Since we used Aeq, the result is ALREADY in the COI frame.
% We double-check just to be safe (removes machine precision noise).
theta_cuep_mod = x_sol - (sum(x_sol .* H_col) / sum(H_col));

% 8. Convergence Check
if exitflag <= 0
    warning('CUEP Optimization did not converge. ExitFlag: %d', exitflag);
    disp(output.message);
else
    [~, mism] = nonlin_con(theta_cuep_mod);
    fprintf('CUEP Found. Max Power Mismatch: %.2e p.u.\n', max(abs(mism)));
    fprintf('Initial Guess Norm: %.4f | Final Solution Norm: %.4f\n', norm(x0), norm(theta_cuep_mod));
end
%%
Vcr_candidate = Calculate_PE_single_point(theta_cuep_mod, ths, Pi, C, D, g);
% ... [Rest of your code] ...
%theta_cuep_mod = [-0.7451 2.3909 0.6487]';

%Vcr_candidate = Calculate_PE_single_point(theta_cuep_mod, ths, Pi, C, D, g);
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