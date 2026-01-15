clc; clear;
%% 1. Setup Data
tic
load('data1.mat'); 
load('Y_all.mat');
load('CCT_TimeDomain.mat', 'CCT_TD');
load('Offline_Database.mat', 'TEF_Database');
load('Tcl.mat', 't_cl');
load('Fault_Info.mat');
 load('shedable_gen.mat');
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
%% COI reference
d_COI = (delta * M) / M_tot;
w_COI = (omega * M) / M_tot;
theta = delta - d_COI;
w_tilde = omega - w_COI;
w = w_tilde;
th=theta;

%% Solving for post fault SEP using Load Flow
idx_gen=[1 2 3 4 5 6 7 8 9 10];
idx_load=setdiff(1:num_bus,idx_gen);
v_gen=[1.0475 0.9820 0.9831 0.9972 1.0123 1.0493 1.0635 1.0278 1.0265 1.0300];
slack_bus=31;
%gen_buses = [30, 31,    32,  33,    34,  35,  36,   37,    38,   39] 
xd_p=     [0.025, 0.05, 0.045, 0.035, 0.089, 0.04, 0.044, 0.045, 0.045, 0.004]';
Pgen=zeros(num_bus,1);
Pgen(1:num_gen)=Pm(1,1:num_gen);
YN_post=Y_post;
x_eq_post=NR_ss(YN_post,Pgen,idx_load,v_gen,slack_bus);
V_eq_post=x_eq_post(num_bus+1:end).*(cos(x_eq_post(1:num_bus))+1i*sin(x_eq_post(1:num_bus)));
I_eq_post=YN_post*V_eq_post;
Pgen_post=real(V_eq_post(idx_gen).*conj(I_eq_post(idx_gen)));
Eeq_post=abs(V_eq_post(idx_gen)+1i*xd_p.*I_eq_post(idx_gen));
delta_eq_post=angle(V_eq_post(idx_gen)+1i*xd_p.*I_eq_post(idx_gen));
x_eq_post=[delta_eq_post-M'*delta_eq_post/M_tot; zeros(num_gen,1)];
%ths=[-0.1782 0.5309 0.2711]';
%% krons reduced matrix
 Y1=Yint_post;
 % E=Eeq_post;
 %E=[1.054 1.050 1.017];
 C=zeros(g,g);
 D=zeros(g,g);
 for i=1:g
    for j=1:g
         C(i,j)=E(i)*E(j)*imag(Y1(i,j));
        D(i,j)=E(i)*E(j)*real(Y1(i,j)); 
    end
 end
Pi = Pgen(1:num_gen) - (real(diag(Y1))) .* ((E(:)).^2);


%% 3. SIMULATE "REAL-TIME" DATA ACQUISITION
fault_start_time = 1.0;
data_time = fault_start_time + CCT_TD + 0.1; % Example: Test near the CCT found earlier

[~, idx_online] = min(abs(T - data_time));

t_online = T(idx_online);
fprintf('2. Acquiring Data Snapshot at t = %.4fs\n', t_online);

% Get Measurements at this instant
delta_online = theta(idx_online, :);
omega_online = w(idx_online, :);

KE=0.5.*M'.*(omega_online).^2;
KE_tot=sum(KE);
KE_Norm=KE/KE_tot;
%% kNN search
X =vertcat([TEF_Database.KE_Signature].');
Y = KE_Norm;
[idx_knn,Euc_D]=knnsearch(X,Y,"K",5,"Distance","euclidean");

Fault = {TEF_Database(idx_knn).Fault_Location}.';
Sorted_Mod={TEF_Database(idx_knn).MOD_Generators};
%% Solving for post fault sep
ths=x_eq_post(1:num_gen);
%% Finding CUEP vector For each ranked MOD
CUEP_Candidate=zeros(length(Sorted_Mod),num_gen);
for k=1:length(Sorted_Mod)
MOD=Sorted_Mod{k};
theta_u = Calculate_theta_u_online(MOD, num_gen, ths, H);
Pm_row = Pm(1, 1:num_gen);   % 1x3 Row Vector
E_row  = E;         % 1x3 Row Vector
H_col  = H(:);               % 3x1 Column Vector
% objective function
obj_fun = @(x) 1; 
%initial condition
x0 = theta_u;
%linear Constraint
A = []; 
b = []; 
Aeq = [];       
beq = []; 

%bounds

lb = -2*pi * ones(g, 1); 
ub =  2*pi * ones(g, 1);


opts = optimset('Algorithm', 'interior-point');
%opts = optimset('Diplay','iter','Algorithm','sqp');
% 5. Constraint Function Handle
nonlin_con = @(x) SEPfunction(x, Pm_row, E_row, C, D, H, Y1);

% 6. Run Optimization
fprintf('Solving for CUEP with COI Constraint...\n');
[x_sol, fval, exitflag, output] = fmincon(obj_fun, x0, A, b, Aeq, beq, lb, ub, nonlin_con, opts);
%%
% 7. Post-Process
theta_cuep_mod= x_sol- (sum(x_sol .* H_col) / sum(H_col));

% 8. Convergence Check
if exitflag <= 0
    warning('CUEP Optimization did not converge. ExitFlag: %d', exitflag);
else
    [~, mism] = nonlin_con(theta_cuep_mod);
    fprintf('CUEP Found. Max Power Mismatch: %.2e p.u.\n', max(abs(mism)));
end
  

fprintf('Saving the controlling Unstable Equillibrium Point \n');
CUEP_Candidate(k,:)=theta_cuep_mod;

end
%% Correct MOD Identification
%Checking for few MODs
num_candidates = length(Sorted_Mod);
limit = min(num_candidates, 2);
Vcr_candidate = zeros(limit, 1);
KE_corr_candidate = zeros(limit, 1);
V_total = zeros(limit, 1);
del_VPE_candidate = zeros(limit, 1);
norm_delVPE = zeros(limit, 1);
VPE= Calculate_PE_single_point(delta_online, ths, Pi, C, D, g);
for i=1:limit
MOD=Sorted_Mod{i};
Vcr_candidate(i) = Calculate_PE_single_point(CUEP_Candidate(i,:), ths, Pi, C, D, g);
KE_corr_candidate(i) = Calculate_KE_online(g, H, Ws, omega_online, MOD);
V_total(i) = VPE + KE_corr_candidate(i);
del_VPE_candidate(i)=Vcr_candidate(i)-V_total(i);
norm_delVPE(i)=del_VPE_candidate(i)/KE_corr_candidate(i);
end

ResultTable = table((1:limit)', Vcr_candidate, KE_corr_candidate, del_VPE_candidate, norm_delVPE, ...
    'VariableNames', {'Rank', 'V_CUEP', 'KE_Corr', 'Margin_DeltaV', 'Norm_Margin'});

disp(ResultTable);
% Find the lowers normalized potential energy
[min_maargin, idx_mod]=min(norm_delVPE);
Correct_MOD=Sorted_Mod{idx_mod};

%% Assessment of the First Swing Stability
 %extracting angle and speed at fault clearing
 faultclearing_time=1+t_cl;
 [~, idx_clearing] = min(abs(T - faultclearing_time));

 wcl=w(idx_clearing, :);
 thetacl=theta(idx_clearing, :);
 thcuep=CUEP_Candidate(idx_mod,:);
 Vcr= Calculate_PE_single_point(thcuep, ths, Pi, C, D, g);
 VPEcl= Calculate_PE_single_point(thetacl, ths, Pi, C, D, g);
 KE_corr_cl = Calculate_KE_online(g, H, Ws, wcl,Correct_MOD);
 Vcl=VPEcl+KE_corr_cl;
 del_V=Vcr-Vcl;
 fprintf("Energy margin=%d\n",del_V);
 
if del_V<0
    Shed=1;
    fprintf("Energy margin is negative System is Unstable take control action \n");
else
    Shed=0;
    fprintf("Energy margin is positive System is Stable no Need to take control Action\n");
end
 %% Shedding Calculation
if Shed==1
 Gen_idx=Gen_sh;
 threshold=0;
 Psh_max=0.40;
 max_iter=10;
data_time=1+t_cl+0.2;
[~, idx_online] = min(abs(T - data_time));
wpred=w(idx_online,:);
deltapred=theta(idx_online,:);

k=1;
VPE(k)=Calculate_PE_single_point(deltapred, ths, Pi, C, D, g);
VKE(k)=Calculate_KE_online(g, H, Ws, w, MOD);
Vtot(k)=VPE(k)+VKE(k);
Vsh(k)=Vcr-Vtot(k);
Pmin=0.05;
k=k+1;
Psh(k)=Pmin;
H1=H;
Pgen1=Pgen;
xd_p1=xd_p;
Hnew=H1(Gen_idx)*(1-Psh(k));
H1(Gen_sh)=Hnew;
Pgen_new=Pgen1(Gen_idx)*(1-Psh(k));
Pgen1(Gen_idx)=Pgen_new;
Xdnew=xd_p1(Gen_idx)/(1-Psh(k));
xd_p1(Gen_idx)=Xdnew;
Pi1 = Pgen1(1:num_gen) - (real(diag(Y1))) .* ((E(:)).^2);
VPE(k)=Calculate_PE_single_point(deltapred, ths, Pi1, C, D, g);
VKE(k)=Calculate_KE_online(g, H1, Ws, w, MOD);
Vtot(k)=VPE(k)+VKE(k);
Vsh(k)=Vcr-Vtot(k);
optimal_shedding_found = false;
final_cmd=0;

while k<max_iter
    if Vsh(k)<threshold
        if Psh(k)<Psh_max
            final_cmd=Psh(k);
            optimal_shedding_found=true;
            break
        else
                warning('Required shedding exceeds limit. Initiating other controls.');
                optimal_shedding_found = false;
                break;
        end

    else
        % "Interpolate/Extrapolate Psh value at V=0"
        % Secant Method for Interpolation:
            % Finds P_new where V approaches 0 based on slope of previous two points
            y2 = Vsh(k);     x2 = Psh(k);
            y1 = Vsh(k-1);   x1 = Psh(k-1);
            slope = (y2 - y1) / (x2 - x1);
            P_new = x2 + (threshold - y2) / slope;
            if P_new <= x2
                P_new = x2 + 0.05; % Force 5% increment if math fails
            end
            k=k+1;
            Psh(k)=P_new;
            H1=H;
            Pgen1=Pgen;
            xd_p1=xd_p;
            Hnew=H1(Gen_idx)*(1-Psh(k));
            H1(Gen_idx)=Hnew;
            Pgen_new=Pgen1(Gen_idx)*(1-Psh(k));
            Pgen1(Gen_sh)=Pgen_new;
            Xdnew=xd_p1(Gen_idx)/(1-Psh(k));
            xd_p1(Gen_idx)=Xdnew;
            Pi1 = Pgen1(1:num_gen) - (real(diag(Y1))) .* ((E(:)).^2);
            VPE(k)=Calculate_PE_single_point(deltapred, ths, Pi1, C, D, g);
            VKE(k)=Calculate_KE_online(g, H1, Ws, w, MOD);
            Vtot(k)=VPE(k)+VKE(k);
            Vsh(k)=Vcr-Vtot(k);
    end
end
if optimal_shedding_found
        fprintf('>>> SUCCESS: Shed %.2f%% of Gen %d to stabilize.\n', final_cmd*100, Gen_sh);
        % Apply actual shedding logic here for simulation...
 else
        fprintf('>>> FAILURE: Cannot stabilize with shedding alone. Take other actions.\n');
 end
end
 Online_comp=toc