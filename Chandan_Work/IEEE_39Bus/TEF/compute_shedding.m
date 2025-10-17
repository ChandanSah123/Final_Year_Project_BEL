% compute_shedding.m
% Given a chosen unstable MOD, compute a candidate generation shedding plan that restores margin.
% Simple approach: use precomputed sensitivity in DB or binary search on P_shed to make Vcl < Vcr.
function controlPlan = compute_shedding(chosenResult, meas, sys)
% chosenResult: a struct returned by TEF_online_assess 'chosen' entry
% For now, return a simple plan: shed 10% from largest generator in MOD
MOD = chosenResult.MOD;
% pick largest Pm generator in MOD
Pm_vals = arrayfun(@(i) sys.gen(i).Pm, MOD);
[~, imax] = max(Pm_vals);
shedGen = MOD(imax);
shedFraction = 0.1;
controlPlan.gen = shedGen;
controlPlan.Pshed = shedFraction * sys.gen(shedGen).Pm;
controlPlan.estimate = 'simple heuristic - replace with sensitivity-based root find';
end
