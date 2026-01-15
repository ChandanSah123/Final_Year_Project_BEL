clc; clear; close all;
% 1. Setup Data
global Yint_post;
load('data1.mat'); 
load('Y_all.mat');
load('CCT_TimeDomain.mat', 'CCT_TD')
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
 plot(fth(1001:1500));

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

%% Finding mod difference (MOD) technique
MOD_idx=2;
MOD=MOD_sort_data(1:MOD_idx,1);
theta_u = Calculate_theta_u(MOD_idx, MOD_sort_data, num_gen, ths, H);


%% Minimization MOD
tic
options = odeset('MaxStep',0.01,'InitialStep',0.01);

% PASS PARAMETERS HERE TOO
[Tm, Y] = ode45(@(t,y) Integrateth(t, y, Pi, C, D, H), [0 9], theta_u, options);

ym = Y;
toc
theta_cuep_mod = ym(end, :)'; 
fprintf('CUEP Angles (radians) MOD: \n');
disp(theta_cuep_mod);

%%
fprintf('Energy margin calculation \n');
 th_vec=theta_cuep_pebs;
 Vcr = Calculate_PE_single_point(th_vec, ths, Pi, C, D, g);
 VPE=Calculate_PE(npts, g, Pi, C, D, th, ths);
 [KE_total, KE_corr] = Calculate_KE(npts, g, H, Ws, w, MOD);
 V=VPE+KE_corr;
 delV=Vcr-V;
 figure; plot(V,'g');
hold on;
 plot(delV,'b');

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

