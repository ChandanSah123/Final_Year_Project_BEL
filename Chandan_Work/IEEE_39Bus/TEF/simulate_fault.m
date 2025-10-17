% simulate_fault.m
% Time-domain simulator wrapper. Uses a classical 2nd-order machine model by default.
% Outputs:
% - traj: struct with full time-series theta(t), omega(t), t
% - meas: struct with theta_cl (angles at clearing), omega_cl, time-series for measurement window
function [traj, meas] = simulate_fault(sys, faultInfo, simOpts)
% simOpts.tfinal, simOpts.dt
tspan = 0:simOpts.dt:simOpts.tfinal;
nGen = numel(sys.gen);
% initial conditions (angles and speeds) from power flow
theta0 = zeros(nGen,1); omega0 = zeros(nGen,1);
% TODO: set initial rotor angles to pre-fault steady-state (solve power flow + swing eq steady)
x0 = [theta0; omega0];

% Define ODE function (classical swing, electrical power computed from network Ybus)
odefun = @(t,x) swing_odes(t,x,sys,faultInfo);
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
[T,X] = ode45(odefun,tspan,x0,options);

theta = X(:,1:nGen);
omega = X(:,nGen+1:end);

traj.t = T;
traj.theta = theta;
traj.omega = omega;

% find clearing index
tclear = faultInfo.t_clear;
[~, idx_clear] = min(abs(T - tclear));
meas.theta_cl = theta(idx_clear,:)';
meas.omega_cl = omega(idx_clear,:)';
% also return measurement time-series for a short window after fault
window = min(find(T >= (faultInfo.t_fault + simOpts.dt),1) : find(T <= (faultInfo.t_fault + 0.12),1,'last'));
if isempty(window)
    window = 1:min(ceil(0.08/simOpts.dt),length(T));
end
meas.t = T(window);
meas.theta_ts = theta(window,:);
meas.omega_ts = omega(window,:);
end

% internal swing ODEs: classical 2nd-order machine model (placeholder)
function dx = swing_odes(t,x,sys,faultInfo)
nGen = numel(sys.gen);
theta = x(1:nGen);
omega = x(nGen+1:end);
% Detect if fault is active: if t between faultInfo.t_fault and t_clear, use faulted Ybus
if t >= faultInfo.t_fault && t <= faultInfo.t_clear
    Y = apply_fault_to_Y(sys.Ybus, faultInfo);
else
    Y = sys.Ybus;
end
% compute electrical powers Pe using internal EMFs and Y (structure-preserving reduction or classical)
Pe = zeros(nGen,1);
for i=1:nGen
    bus_i = sys.gen(i).bus;
    Ei = sys.gen(i).E;
    % approximate: Pe = sum_j |Ei||Ej| B_ij sin(theta_i - theta_j)
    % TODO: convert to more accurate network calculation using internal EMFs and reduced Ybus
    for j=1:nGen
        if i==j, continue; end
        Ej = sys.gen(j).E;
        Bij = imag(-1/sys.Ybus(bus_i,bus_i)); % placeholder (WRONG for real networks) -> replace
        Pe(i) = Pe(i) + Ei*Ej*0.0; % placeholder zero; replace with correct network computation
    end
    Pe(i) = sys.gen(i).Pm; % placeholder to avoid dynamics until replaced
end

% Swing equations: M * d(omega) = Pm - Pe - D*omega
dtheta = omega;
domega = zeros(nGen,1);
for i=1:nGen
    M = sys.gen(i).M;
    D = sys.gen(i).D;
    domega(i) = (sys.gen(i).Pm - Pe(i) - D*omega(i))/M;
end
dx = [dtheta; domega];
end

% NOTE: simulate_fault uses placeholders for electrical power calculation.
% Replace apply_fault_to_Y and Pe calc with correct reduced-order classical model or external simulator.
