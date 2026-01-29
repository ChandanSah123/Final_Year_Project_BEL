clc; clear; close all;
tic
load('data1.mat'); 
load('Y_all.mat');
load('CCT_TimeDomain.mat', 'CCT_TD');
load('Tcl.mat')
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
E= [1.0562; 1.1782; 1.1507; 1.1004; 1.3594; 1.1806; 1.1288; 1.0721; 1.1382; 1.0215]';
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
idx_gen=[1 2 3 4 5 6 7 8 9 10];
idx_load=setdiff(1:num_bus,idx_gen);
Pgen=zeros(num_gen,1);
Pgen(1:num_gen)=Pm(1,1:num_gen);
%ths=x_eq_post(1:num_gen);
%ths=[-0.1782 0.5309 0.2711]';

%% krons reduced matrix
  Y1=Yint_post;
  %E=Eeq_post;
 C=zeros(g,g);
 D=zeros(g,g);
 for i=1:g
    for j=1:g
         C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j)); 
    end
 end
%% Alternative method to find the stable equilibrium point
Pm = Pm(1, 1:num_gen);   % 1x3 Row Vector
%E = Eeq_post.';         % 1x3 Row Vector
H = H(:); 
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


opts = optimset('Algorithm', 'sqp');
%opts = optimset('Diplay','iter','Algorithm','sqp');
% 5. Constraint Function Handle
nonlin_con =@(x) SEPfunction(x, Pm, E, C, D, H, Y1);

% 6. Run Optimization
fprintf('Solving for CUEP with COI Constraint...\n');
[x_sol1, fval1, exitflag1, output1] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, nonlin_con, opts);
ths = x_sol1-(sum(x_sol1(:) .* H(:)) / sum(H));
%display(ths);  -0.183008274905452	0.544513400978037	0.279544801496823
%display(ths1);
 

%% 3. DIRECT CUEP & MOD IDENTIFICATION (Using CCT Snapshot)
fault_start_time = 1.0; 
t_cct_absolute1 = fault_start_time + CCT_TD+0.2; %slightly higher than cct
t_cct_absolute=fault_start_time + CCT_TD;

fprintf('Fault Duration: %.4fs | Absolute Clearing Time: %.4fs\n', CCT_TD, t_cct_absolute1);

[~, idx_cct] = min(abs(T - t_cct_absolute));

fprintf('Snapshot taken at simulation step: %d (t = %.4fs)\n', idx_cct, T(idx_cct));

theta_guess = theta(idx_cct, :)';
[sorted_angles, sort_idx] = sort(theta_guess, 'descend');
w_guess = w_tilde(idx_cct, :)'; 
sorted_speeds = w_guess(sort_idx);
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
tic
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
Aeq = H';       
beq = 0; 

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
%%
% 7. Post-Process
theta_cuep_mod = x_sol - (sum(x_sol(:) .* H_col(:)) / sum(H_col));

% 8. Convergence Check
if exitflag <= 0
    warning('CUEP Optimization did not converge. ExitFlag: %d', exitflag);
else
    %[~, mism] = nonlin_con(theta_cuep_mod);
   % fprintf('CUEP Found. Max Power Mismatch: %.2e p.u.\n', max(abs(mism)));
end
% ... [Rest of your code] ...
%theta_cuep_mod = [-0.7451 2.3909 0.6487]';

Vcr_candidate = Calculate_PE_single_point(theta_cuep_mod, ths, Pi, C, D, g);
%Vcr_candidate = Calculate_PE(1, g, Pi, C, D, theta_u, ths);
[KE, KE_corr_candidate] = Calculate_KE(npts, g, H, Ws, w, current_MOD_indices);
V_total = VPE + KE_corr_candidate;
delV_candidate = Vcr_candidate - V_total;
%plot(delV_candidate);
unstable_indices = find(delV_candidate < 0);

if isempty(unstable_indices)
    fprintf('System is stable for the entire duration of T.\n');
    tcr=t_cl;
    
else
    idx_unstable = unstable_indices(1);
    if idx_unstable == 1
        tcr = 0;
        warning('System is unstable from the start.');
    else
        t_stable   = T(idx_unstable - 1);
        t_unstable = T(idx_unstable);
        margin_stable   = delV_candidate(idx_unstable - 1);
        margin_unstable = delV_candidate(idx_unstable);
        % Linear Interpolation formula to find where Margin = 0
        slope = (margin_unstable - margin_stable) / (t_unstable - t_stable);
        dt = (0 - margin_stable) / slope;
        tcr = t_stable + dt;
        fprintf('CCT Found: %.4fs (Interpolated)\n', tcr);
    end
end
display(tcr);
toc
display(ths);
display(theta_u);
display(theta_cuep_mod);
fprintf("Corrected Kinetic Energy=\n");
display(KE_corr_candidate(idx_cct));
fprintf("Potential Energy=\n");
display(VPE(idx_cct));
fprintf("Total Energy=\n");
display(V_total(idx_cct));
display(Vcr_candidate);
fprintf("Energy margin=\n");
display(delV_candidate(idx_cct));

%% calculation of ratio of kinetic energy to power


%% 4. COMPREHENSIVE ENERGY PLOTTING

figure('Name', 'TEF Analysis: Energy & Stability Margin', 'Color', 'w', 'Position', [100, 100, 1000, 800]);

% --- Subplot 1: System Energies ---
subplot(2,1,1);
hold on;

% Plot Energy Components
plot(T, VPE, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Potential Energy (V_{PE})');
plot(T, KE_corr_candidate, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Corrected KE (V_{KE|corr})');
plot(T, V_total, 'k-', 'LineWidth', 2, 'DisplayName', 'Total Energy (V_{total})');

% Plot Critical Energy Threshold (Horizontal Line)
yline(Vcr_candidate, 'r--', 'LineWidth', 2, 'DisplayName', ['Critical Energy (V_{cr}=' num2str(Vcr_candidate, '%.4f') ')']);
% Mark the CCT (Vertical Line)
xline(tcr, 'm-.', 'LineWidth', 1.5, 'DisplayName', ['CCT = ' num2str(tcr, '%.4f') 's']);

% Aesthetics
ylabel('Energy (p.u.)', 'FontSize', 12, 'FontWeight', 'bold');
title('Transient Energy Components vs. Critical Threshold', 'FontSize', 14);
legend('Location', 'bestoutside');
grid on;
grid minor;
box on;
xlim([min(T) max(T)]); % Adjust x-limit if you want to zoom in

% --- Subplot 2: Energy Margin (Delta V) ---
subplot(2,1,2);
hold on;

% Plot Energy Margin
% Using 'area' makes the positive (stable) and negative (unstable) regions clearer
area(T, delV_candidate, 'FaceColor', [0.85 0.9 1], 'EdgeColor', 'b', 'LineWidth', 1.5, 'DisplayName', 'Energy Margin (\Delta V)');

% Plot Zero Line (Stability Boundary)
yline(0, 'k-', 'LineWidth', 2, 'DisplayName', 'Stability Boundary');

% Mark the CCT point explicitly
plot(tcr, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Instability Point');
xline(tcr, 'm-.', 'LineWidth', 1.5, 'HandleVisibility', 'off'); % Vertical line reference

% Aesthetics
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy Margin (p.u.)', 'FontSize', 12, 'FontWeight', 'bold');
title(['Stability Margin (\Delta V = V_{cr} - V_{total}) | Crossing at t = ' num2str(tcr, '%.4f') 's'], 'FontSize', 14);
legend('Location', 'bestoutside');
grid on;
grid minor;
box on;
xlim([min(T) max(T)]);

% Text Annotations for Stable/Unstable Regions
yl = ylim; 
text(min(T) + 0.02, yl(2)*0.8, 'STABLE REGION (\Delta V > 0)', 'Color', [0 0.5 0], 'FontWeight', 'bold');
text(min(T) + 0.02, yl(1)*0.8, 'UNSTABLE REGION (\Delta V < 0)', 'Color', [0.8 0 0], 'FontWeight', 'bold');

hold off;
