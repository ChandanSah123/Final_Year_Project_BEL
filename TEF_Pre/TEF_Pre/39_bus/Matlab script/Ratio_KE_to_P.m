clc;clear;
%% 
load('data1.mat');
load('Tcl.mat');
num_bus = 39;
num_gen = 10;
g = num_gen;
H = [42.0;      30.3;   35.8;     28.6;   26.0;     34.8;  26.4;       24.3;   34.5;       500.0];
M = 2 * H / (2 * pi * 60); 
M_tot = sum(M);
Ws = (2 * pi * 60);
T = data(:,51);
idx_ang = 1:10;       % Columns for Angle
idx_spd = 11:20;      % Columns for Speed
idx_pm  = 21:30;      % Columns for Mech Power
idx_pe  = 31:40;
delta = data(:, idx_ang) * (pi/180);
omega = data(:, idx_spd) * Ws; 
Pm    = data(:, idx_pm);
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
