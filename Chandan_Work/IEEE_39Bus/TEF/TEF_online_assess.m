% TEF_online_assess.m
% Main online routine implementing the Bhui & Senroy flow using DB.
function [stable, chosenMOD, Vcl, Vcr, controlPlan] = TEF_online_assess(meas, faultInfo, DB, sys, opts)
% meas: struct from simulate_fault (or PMU stream)
% faultInfo: fault definition
% DB: offline lookup database array
% opts.window: measurement window (s)
% opts.kNN_K: number of k-NN neighbors to test
% Returns:
% stable: boolean
% chosenMOD: vector of generator indices forming MOD
% Vcl, Vcr: numeric energy values
% controlPlan: struct with suggested control action (if unstable)
if nargin<5, opts = struct('window',0.08,'kNN_K',5); end

% 1) update post-fault Ybus (fast update)
Ypost = util_Ybus_update(sys.Ybus, faultInfo);
topHash = topology_hash(Ypost);

% 2) compute measured KE_norm from meas window
KE = compute_ke_from_meas(meas, sys, opts.window);
KE_norm = KE / (sum(KE)+eps);

% 3) rank candidates from DB
cand_idxs = rank_MODs_kNN(KE_norm, DB, faultInfo.location, topHash, opts.kNN_K);
if isempty(cand_idxs)
    error('No DB entries matching this fault location/topology');
end

results = struct();
for k = 1:numel(cand_idxs)
    entry = DB(cand_idxs(k));
    % compute post-fault SEP (cheap Newton)
    theta_s = compute_postfault_SEP(Ypost, sys);
    % get CUEP: if DB entry topology matches exactly, reuse theta_u else compute
    if strcmp(entry.topologyHash, topHash)
        theta_u = entry.theta_u;
        Vcr_k = entry.Vcr;
    else
        theta_u = compute_CUEP_via_corner_and_refine(entry.MOD, theta_s, sys, Ypost);
        Vcr_k = compute_energy(theta_u, zeros(size(theta_u)), theta_s, sys);
    end
    % compute Vcl at clearing instant using meas.theta_cl and meas.omega_cl
    theta_cl = meas.theta_cl; omega_cl = meas.omega_cl;
    Vcl_k = compute_energy(theta_cl, omega_cl, theta_s, sys);
    % compute Bhui metric (potential difference normalized by corrected KE)
    KEcorr = sum(KE); % placeholder: implement two-machine equivalent corrected KE
    VPE_u = compute_VPE(theta_u, theta_s, sys);
    VPE_cl = compute_VPE(theta_cl, theta_s, sys);
    metric = (VPE_u - VPE_cl) / (KEcorr + eps);
    results(k).DBidx = cand_idxs(k);
    results(k).MOD = entry.MOD;
    results(k).Vcr = Vcr_k;
    results(k).Vcl = Vcl_k;
    results(k).metric = metric;
end

% choose candidate with minimum metric
metrics = [results.metric];
[~, iminx] = min(metrics);
chosen = results(iminx);
chosenMOD = chosen.MOD;
Vcl = chosen.Vcl;
Vcr = chosen.Vcr;
stable = (Vcl < Vcr);

controlPlan = [];
if ~stable
    % compute suggested generation shedding using precomputed sensitivities if present
    controlPlan = compute_shedding(chosen, meas, sys);
end
end
