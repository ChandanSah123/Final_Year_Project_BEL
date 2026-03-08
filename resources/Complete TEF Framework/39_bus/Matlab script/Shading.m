clc;clear;
%% load the predicted generated angle and speed  from lstm say,
load('data1.mat');
load('Tcl.mat');
load('shedable_gen.mat');
num_bus = 9;
num_gen = 3;
g = num_gen;
H = [23.64; 6.4; 3.01]; 
M = 2 * H / (2 * pi * 60); 
xd_p=[0.0608;0.1198;0.1813];
M_tot = sum(M);
Ws = (2 * pi * 60);
T = data(:,16);
idx_ang = 1:3;
idx_spd = 4:6;
idx_pm  = 7:9;
idx_pe  = 10:12;
delta = data(:, idx_ang) * (pi/180);
omega = data(:, idx_spd) * Ws; 
Pm    = data(:, idx_pm);
%% COI reference
d_COI = (delta * M) / M_tot;
w_COI = (omega * M) / M_tot;
theta = delta - d_COI;
w_tilde = omega - w_COI;
w = w_tilde;
th=theta;
% lets supposet we have to predict the rotor angle at the time instant 0.2
% second after fault clearing
data_time=1+t_cl+0.2;
[~, idx_online] = min(abs(T - data_time));
wpred=w(idx_online,:);
deltapred=theta(idx_online,:);

k=1;
VPE=Calculate_PE_single_point(deltapred, ths, Pi, C, D, g);
VKE=Calculate_KE_online(g, H, Ws, w, MOD);
V=VPE+VKE;
Vsh(k)=V;
Psh(k)=0;
k=k+1;
% lets shed 20 percentage as the minimum shedding from plantsh
shed_p=0.2;

