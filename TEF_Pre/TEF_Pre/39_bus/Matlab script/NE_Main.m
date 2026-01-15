clc; clear; close all;
% 1. Setup Data
global Yint_post;
load('data1.mat'); 
load('Y_all.mat')

T = data(:,51);
npts=length(T);
num_bus=39;
num_gen=10;
g=num_gen;
idx_ang = 1:10;       % Columns for Angle
idx_spd = 11:20;      % Columns for Speed
idx_pm  = 21:30;      % Columns for Mech Power
idx_pe  = 31:40;      % Columns for Elec Power

% Parameters
H = [42.0;      30.3;   35.8;     28.6;   26.0;     34.8;  26.4;       24.3;   34.5;       500.0]; % Example
M = 2 * H / (2 * pi * 60); 
M_tot = sum(M);
Ws=(2*pi*60);

% ... (After loading data)

g = num_gen;
% Now extract data (It will now be perfectly aligned with H and Y1)
delta = data(:, idx_ang)*(pi/180);
omega = data(:, idx_spd)*Ws; 
Pm    = data(:, idx_pm);
Pe    = data(:, idx_pe);
% -------------------------------------------------------------
%% 2. Efficient COI Calculation (Vectorized)
% Calculate Center of Inertia (COI) without loops
d_COI = (delta * M) / M_tot; % Matrix multiplication (N x 10) * (10 x 1) -> (N x 1)
w_COI = (omega * M) / M_tot;

% Transform to COI frame
theta = delta - d_COI;
w_tilde = omega - w_COI;
w = w_tilde;
%% Solving for post fault SEP using Load Flow
idx_gen = 30:39;
idx_load=setdiff(1:num_bus,idx_gen);
v_gen=[1.0475 0.9820 0.9831 0.9972 1.0123 1.0493 1.0635 1.0278 1.0265 1.0300];
slack_bus=31;
%gen_buses = [30, 31,    32,  33,    34,  35,  36,   37,    38,   39] 
xd_p=     [0.025, 0.05, 0.045, 0.035, 0.089, 0.04, 0.044, 0.045, 0.045, 0.004]';
Pgen=zeros(39,1);
Pgen(30:39)=Pm(1,1:10);
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
 %E=[ 1.0562;    1.1782; 1.1507; 1.1004; 1.3594; 1.1806; 1.1288;  1.0721;   1.1382;   1.0215];
 C=zeros(g,g);
 D=zeros(g,g);
 for i=1:g
    for j=1:g
         C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j)); 
    end
 end
%% Fth formulation
ths = x_eq_post(1:10);
th = theta;

% Calculate Pi using the NOW ALIGNED Pm vector
% Pm is (N x 10), we take the first row (t=0) or mean
Pm_aligned = Pm(1, :)'; % 10x1 vector
Pi = Pm_aligned - (real(diag(Y1))) .* ((E(:)).^2);

fth = Calculate_fth(ths, th, Pi, npts, g, C, D, H);
figure; plot(fth(1000:4000)); title('fth Curve '); grid on;

%% Detect PEBS Crossing (Modified)
% 1. DEFINE START POINT
fault_start_time = 1.0;  % Change this to match your exact fault inception time
t_sim=T;
[MOD_sort_data, t_pebs_idx, PEBS_Crossing] = detect_pebs_crossing(fth, t_sim, fault_start_time, npts, th, w_tilde);
%% Finding mod difference (MOD) technique
MOD_idx=5;
MOD=MOD_sort_data(1:MOD_idx,1);
theta_u = Calculate_theta_u(MOD_idx, MOD_sort_data, num_gen, ths, H);
%% minimizing using power flow
%theta_pre1=0.01*rand(num_bus,1); theta_pre1(slack_bus)=0;
%theta_pre1(30:39)=theta_u;
%x_eq_post1=NR_ss(YN_post,Pgen,idx_load,v_gen,slack_bus,theta_pre1);
%V_eq_post1=x_eq_post1(num_bus+1:end).*(cos(x_eq_post1(1:num_bus))+1i*sin(x_eq_post1(1:num_bus)));
%I_eq_post1=YN_post*V_eq_post1;
%Pgen_post=real(V_eq_post1(idx_gen).*conj(I_eq_post1(idx_gen)));
%Eeq_post=abs(V_eq_post1(idx_gen)+1i*xd_p.*I_eq_post1(idx_gen));
%delta_eq_post1=angle(V_eq_post1(idx_gen)+1i*xd_p.*I_eq_post1(idx_gen));
%x_eq_post1=[delta_eq_post1-M'*delta_eq_post/M_tot; zeros(num_gen,1)];
%theta_cuep_mod=x_eq_post1(idx_ang);
%disp(theta_cuep_mod');

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
disp(theta_cuep_pebs');



%% Minimization MOD
tic
options = odeset('MaxStep',0.01,'InitialStep',0.01);

% PASS PARAMETERS HERE TOO
[Tm, Y] = ode45(@(t,y) Integrateth(t, y, Pi, C, D, H), [0 9], theta_u, options);

ym = Y;
toc
theta_cuep_mod = ym(end, :)'; 
fprintf('CUEP Angles (radians) MOD: \n');
disp(theta_cuep_mod');

%% 1. Robust CUEP Solver (Replaces ode45 for PEBS)
fprintf('Finding CUEP (PEBS) using Constrained Solver...\n');

% Start Point: The crossing point found by PEBS
th_start = PEBS_Crossing(:, 2); 

% Define Bounds: Constrain search to +/- 2.5 radians (~140 deg) from exit
% This PREVENTS the "Pole Slip" result (angles jumping 6.28 rad)
lb = th_start - 2.5; 
ub = th_start + 2.5;

% Objective: Minimize sum of squares of Accelerating Power (Mismatch)
% We reuse your existing logic: Pi - P_elec
calc_mismatch_sq = @(th) (Pi - ( ...
    sum( C .* sin(th * ones(1,g) - ones(g,1) * th') + ...
         D .* cos(th * ones(1,g) - ones(g,1) * th'), 2 ) ...
)).^2;

options = optimoptions('lsqnonlin', 'Display', 'off', ...
    'FunctionTolerance', 1e-6, 'StepTolerance', 1e-6);

% Solve!
[theta_cuep_pebs, resnorm] = lsqnonlin(@(x) sqrt(calc_mismatch_sq(x)), th_start, lb, ub, options);

fprintf('CUEP Angles (radians) PEBS (Fixed): \n');
disp(theta_cuep_pebs');

%% 2. Robust CUEP Solver (Replaces ode45 for MOD)
%fprintf('Finding CUEP (MOD) using Constrained Solver...\n');

% Recalculate theta_u based on sorted data (Your existing function)
%theta_u = Calculate_theta_u(MOD_idx, MOD_sort_data, num_gen, ths, H);

% Bounds for MOD search
%lb_mod = theta_u - 2.5;
%ub_mod = theta_u + 2.5;

% Solve!
%[theta_cuep_mod, resnorm_mod] = lsqnonlin(@(x) sqrt(calc_mismatch_sq(x)), theta_u, lb_mod, ub_mod, options);

%fprintf('CUEP Angles (radians) MOD (Fixed): \n');
%disp(theta_cuep_mod');

%%
fprintf('Calculating Critical energy at unstable equillibrium point \n');
 th_vec=theta_cuep_pebs;
 Vcr = Calculate_PE_single_point(th_vec, ths, Pi, C, D, g);
 VPE=Calculate_PE(npts, g, Pi, C, D, th, ths);
 [KE_total, KE_corr] = Calculate_KE(npts, g, H, Ws, w, MOD);
 V=VPE+KE_corr;
 delV=Vcr-V;
 figure;
 plot(V,'g'); 
hold on;
plot(delV,'b');
hold on;
plot(VPE,'r');
title('V green and DelV blue VPE red');
 

% Corrected Loop to find Critical Clearing Time (CCT)
fault_start_time = 1.0;
t_bcu_idx = 0;
    
% 1. DEFINE START POINT (Ensure we start after fault inception)
idx_start = find(T >= fault_start_time, 1);
if isempty(idx_start), idx_start = 2; end

t = idx_start + 1;
found_crossing = false;

while t <= npts
    % CHECK: Did Margin drop from Positive (Safe) to Negative (Unstable)?
    if delV(t) < 0 && delV(t-1) >= 0
        t_bcu_idx = t;
        found_crossing = true;
        
        % Linear Interpolation for exact time
        y1 = delV(t-1);
        y2 = delV(t);
        fraction = y1 / (y1 - y2); % Distance from y1 to 0
        cct_val = T(t-1) + fraction * (T(t) - T(t-1));
        
        fprintf('CCT determined by Energy Crossing at t = %.4f s (Index %d)\n', cct_val, t_bcu_idx);
        break;
    end
    t = t + 1;
end
