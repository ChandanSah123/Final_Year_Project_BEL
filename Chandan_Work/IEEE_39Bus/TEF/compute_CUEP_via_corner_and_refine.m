% compute_CUEP_via_corner_and_refine.m
% Use Bhui corner-point initial guess and local refinement (Newton / optimization)
function theta_u = compute_CUEP_via_corner_and_refine(MOD, theta_s, sys, Ypost)
n = numel(sys.gen);
% Simple corner-point guess:
theta_u = theta_s;
% Choose separation magnitude delta (initial)
delta = 0.5; % rad initial guess - tune or compute via Bhui eqns (3)-(5)
I = MOD; % machines in unstable group
groupI = false(n,1); groupI(I) = true;
theta_u(groupI) = theta_s(groupI) + delta;
theta_u(~groupI) = theta_s(~groupI) - delta*(sum(groupI)/sum(~groupI));

% Local refinement: minimize norm of electrical mismatch f_i(theta) for i in boundary set
% Use fsolve or trust-region to solve f(theta) = 0 for unstable equilibrium (set omega=0)
opts = optimoptions('fsolve','Display','off','TolFun',1e-8,'TolX',1e-8);
try
    theta_u = fsolve(@(th) UEQ_mismatch(th, theta_s, sys, Ypost, MOD), theta_u, opts);
catch
    warning('fsolve failed to refine CUEP; returning corner-point guess.');
end
end

function F = UEQ_mismatch(theta, theta_s, sys, Ypost, MOD)
% Compute mismatches f_i(theta) = mechanical torque - electrical torque (set omega=0)
n = numel(sys.gen);
F = zeros(n,1);
% TODO: compute per-machine electrical power f_i using network and EMFs
for i=1:n
    Pm = sys.gen(i).Pm;
    Pelec = 0; % compute from Ypost and EMFs
    F(i) = Pm - Pelec;
end
% Optionally only return residuals for machines in the boundary set (others fixed)
end
