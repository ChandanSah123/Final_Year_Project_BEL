% example_main.m
% Demonstration / entry point for IEEE-9 TEF framework skeleton
% Usage: run example_main after filling in system data loader and any TODOs.

clear; clc; close all;

% Add project paths (if functions live in subfolders)
% addpath('src');

%% 1) Load system data (buses, lines, generator params)
sys = read_system_data('ieee9');  % implement loader to return struct with fields below

% sys should include:
% sys.Ybus, sys.bus, sys.line, sys.gen (with M, D, E, Xd', Pm, H etc.), sys.baseMVA

%% 2) Offline DB build (optional heavy step)
opts.offline_build = true;
opts.operating_points = {'nominal'}; % extend as needed
opts.fault_list = {'line_4_5', 'bus_6_fault'}; % list of faults to precompute
opts.offline_dt = 0.1;  % seconds, KE fingerprint window
if opts.offline_build
    DB = TEF_offline_build(sys, opts);
    save('TEF_DB_ieee9.mat','DB');
else
    load('TEF_DB_ieee9.mat','DB');
end

%% 3) Simulate an online fault (for testing) or accept PMU stream
% simulate a fault (example): line_4_5, fault-on at t=0, cleared at t=0.08s
faultInfo.location = 'line_4_5';
faultInfo.t_fault = 0.0;
faultInfo.t_clear = 0.08;

% run time-domain sim to generate 'real' measured data
simOpts.tfinal = 5.0;
simOpts.dt = 1e-3;
[traj, meas] = simulate_fault(sys, faultInfo, simOpts);
% meas should contain theta(t), omega(t) for generator buses and clearing instant data

%% 4) Online TEF assessment
onlineOpts.window = 0.08; % seconds of measurements after onset/clearing used to compute KE
onlineOpts.kNN_K = 5;     % k-NN neighbors to test
[stable, chosenMOD, Vcl, Vcr, controlPlan] = TEF_online_assess(meas, faultInfo, DB, sys, onlineOpts);

fprintf('TEF decision: stable=%d, chosenMOD=%s, Vcl=%.4f, Vcr=%.4f\n', ...
    stable, mat2str(chosenMOD), Vcl, Vcr);
if ~isempty(controlPlan)
    disp('Control plan (example):');
    disp(controlPlan);
end
