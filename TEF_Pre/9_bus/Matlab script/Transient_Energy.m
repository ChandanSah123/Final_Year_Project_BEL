clc; clear; close all;
% 1. Setup Data
global Yint_post;
load('data1.mat'); 
load('Y_all.mat')
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
% E=Eeq_post;
 E=[1.054 1.050 1.017];
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
 fth=zeros(npts,1);
 Pei=zeros(1,g);
 f=zeros(1,g);
 for t=1:npts
    Pcoi=0;
    for i=1:g
        Pcoi=Pcoi+Pi(i);
        for j=i+1:g
            Pcoi=Pcoi-2*D(i,j)*cos(th(t,i)-th(t,j));
        end
    end
    for i=1:g
        Pei(i)=0;
        for j=1:g
            if j~=i
                Pei(i)=Pei(i)+C(i,j)*sin(th(t,i)-th(t,j))+D(i,j)*cos(th(t,i)-th(t,j));
            end
        end
        f(i)=Pi(i)-Pei(i)-(H(i)/sum(H))*Pcoi;
        fth(t)=fth(t)+((th(t,i)-ths(i))*f(i));
    end
 end
% figure;plot(fth,'g');
% plot(fth(1001:1900));

%% Detect PEBS Crossing (Modified)

% 1. DEFINE START POINT
% We must ignore the pre-fault period (e.g., 0s to 1s) to avoid false triggers.
fault_start_time = 1.0;  % Change this to match your exact fault inception time
t_sim=T;
idx_start = find(t_sim >= fault_start_time, 1);

if isempty(idx_start)
    idx_start = 2; % Fallback: start at index 2 if time vector not found
else
    idx_start = idx_start + 1; % Start 1 step after fault to allow comparison
end

t = idx_start + 1;
t_pebs_idx = 0;

% 2. Find the Crossing (Loop through Fault Trajectory only)
while t < npts
    % CRITICAL CHANGE: 
    % We look for a sign change (Negative -> Positive).
    % This filters out the initial positive noise at equilibrium.
    if fth(t) > 0 && fth(t-1) <= 0
        t_pebs_idx = t;
        break;
    end
    t = t + 1;
end

if t_pebs_idx > 0
    % 3. Extract Snapshot (Using Linear Interpolation for better accuracy)
    % Calculate fraction where crossing happened between t-1 and t
    y1 = fth(t_pebs_idx-1);
    y2 = fth(t_pebs_idx);
    fraction = abs(y1) / (abs(y1) + abs(y2));
    
    % Interpolate Angles and Speeds
    raw_angles_prev = th(t_pebs_idx-1, :);
    raw_angles_curr = th(t_pebs_idx, :);
    raw_angles = raw_angles_prev + fraction * (raw_angles_curr - raw_angles_prev);
    
    raw_speeds_prev = w_tilde(t_pebs_idx-1, :);
    raw_speeds_curr = w_tilde(t_pebs_idx, :);
    raw_speeds = raw_speeds_prev + fraction * (raw_speeds_curr - raw_speeds_prev);
    
    gen_ids = 1:length(raw_angles);   
    
    % 4. Sort by Angle (Descending)
    [sorted_angles, sort_idx] = sort(raw_angles, 'descend');
    
    % 5. Match IDs and Speed to the sorted order
    sorted_speeds  = raw_speeds(sort_idx);
    sorted_ids     = gen_ids(sort_idx);
    
    % 6. Create ONE Matrix
    MOD_sort_data = [sorted_ids', sorted_angles', sorted_speeds'];
    
    % 7. Save ONLY this matrix to the file
    save('C_MOD','MOD_sort_data');
    
    fprintf('Success! PEBS Crossing found at t = %.4f s (Index %d).\n', t_sim(t_pebs_idx), t_pebs_idx);
    disp('Variable "MOD_sort_data" content:');
    disp(MOD_sort_data);
else
    disp('No PEBS crossing found. The fault might be cleared too fast or simulation time is too short.');
end%% Determining the mod


%% Finding mod difference (MOD) technique
MOD_idx=2;
MOD=MOD_sort_data(1:MOD_idx,1);
%[dt1, dt2]=MOD_SHIFT(MOD, num_gen, ths,H);
%thu1=ths+dt1;
%thu2=ths-dt2;
group1=MOD;
all_gen=1:num_gen;
group2 = setdiff(all_gen, group1);
M1=sum(M(group1));
M2=sum(M(group2));
MT=M1+M2;
theta1_s=sum(ths(group1).*M(group1))/M1;
theta2_s=sum(ths(group2).*M(group2))/M2;
d_theta_s=theta1_s-theta2_s;
d_theta1=(pi-(2*d_theta_s))*(M2/MT);
d_theta2=(pi-(2*d_theta_s))*(M1/MT);
 
theta_u(group1)=ths(group1)+d_theta1;     % belong to group1
theta_u(group2)=ths(group2)-d_theta2;     % belong to group2

%% performing the gradinet decent/sigma F(th) minimization using pebs data
tic
thpebs= [-0.719547621461578	1.89292330246141	1.62637761980022];
%-0.719547621461578	1.89292330246141	1.62637761980022
options = odeset('MaxStep',0.01,'InitialStep',0.01);
[Tp,Y]=ode45(@Integrateth,[0 3],thpebs,options);
yp=Y;
toc
theta_cuep_pebs = yp(end, :)';  % Transpose to make it a column vector (3x1)

fprintf('CUEP Angles (radians) BCU: \n');
disp(theta_cuep_pebs);
%% performing the gradinet decent/sigma F(th) minimization using pebs data
tic

options = odeset('MaxStep',0.01,'InitialStep',0.01);
[Tm,Y]=ode45(@Integrateth,[0 3],theta_u,options);
ym=Y;
toc
theta_cuep_mod = ym(end, :)';  % Transpose to make it a column vector (3x1)

fprintf('CUEP Angles (radians) MOD: \n');

disp(theta_cuep_mod);

%%
  KE=zeros(npts,1);
  KE1=zeros(npts,1);
  PE1=zeros(npts,1);
  PE2=zeros(npts,1);
  PE3=zeros(npts,1);
  PE=zeros(npts,1);
  V=zeros(npts,1);

  wc=w_tilde;
  for t=1:npts
    for i=1:g
        KE1(t)=KE1(t)+0.5*(2*H(i)/Ws)*(wc(t,i))^2;
        PE1(t)=PE1(t)+Pi(i)*(th(t,i)-ths(i));
    end
    wcr=0;
    wsys=0;
    Hcr=0;
    Hsys=0;
   % MOD=MOD1;
    for i=1:g
    c=0;
    for j=1:numel(MOD)
        if i==MOD(j)
            wcr=wcr+H(i)*wc(t,i);
            Hcr=Hcr+H(i);
            c=1;
        end
     end
    if c==0
        wsys=wsys+H(i)*wc(t,i);
        Hsys=Hsys+H(i);
    end
    end
    wcr=wcr/Hcr;
    wsys=wsys/Hsys;
    weq(t)=wcr-wsys;
    Heq=Hcr*Hsys/(Hcr+Hsys);
    KE(t)=0.5*(2*Heq/Ws)*(weq(t))^2;
    for i=1:g-1
        for j=i+1:g
            PE2(t)=PE2(t)+C(i,j)*(cos(th(t,i)-th(t,j))-cos(ths(i)-ths(j)));
            PE3(t)=PE3(t)-D(i,j)*((th(t,i)-ths(i)+th(t,j)-ths(j))/(th(t,i)-th(t,j)-ths(i)+ths(j)))*(sin(th(t,i)-th(t,j))-sin(ths(i)-ths(j)));
        end
    end
    PE(t)=(PE1(t)+PE2(t)+PE3(t));
    KE(t)=KE(t);
    V(t)=KE(t)+PE(t); %for PEBS detection KE-PE, for calculating dV.. -KE+PE
  end
  plot(V);