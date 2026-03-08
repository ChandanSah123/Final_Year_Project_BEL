clc;clear;
%% 
load('data1.mat');
load('Tcl.mat');
T = data(:,26);
npts=length(T);
num_bus=14;
num_gen=5;
g=num_gen;
idx_ang = 1:5;       % Columns for Angle
idx_spd = 6:10;      % Columns for Speed
idx_pm  = 11:15;      % Columns for Mech Power
idx_pe  = 16:20;      % Columns for Elec Power

% Parameters
H = [5.148; 6.54; 6.54; 5.06; 5.06]; % Example
E=[1.1188 1.0893 1.0418 1.0706 1.0970];
xd_p = [0.230; 0.13; 0.13; 0.12; 0.12];
M = 2 * H / (2 * pi * 60); 
M_tot = sum(M);
Ws=(2*pi*60);
% Extract Data
delta = data(:, idx_ang)*(pi/180);
omega = data(:, idx_spd)*Ws;
Pm = data(:, idx_pm);
Pe = data(:, idx_pe);
% accessing data at clearing time
data_time=1+t_cl;
[~, idx_online] = min(abs(T - data_time));
wpred=omega(idx_online,:);
deltapred=delta(idx_online,:);
Pgen=Pm(100,1:num_gen);
KE=0.5*M'.*(wpred).^2;
Ratio= KE ./ Pgen;
[max_val, idx] = max(Ratio);
Gen_sh=idx;
save('shedable_gen.mat','Gen_sh','Pgen');
