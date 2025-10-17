% compute_ke_from_meas.m
% Compute per-generator kinetic energy given measurement time-series meas (theta_ts, omega_ts)
% Returns vector ke(i) computed over the window (use difference from pre-fault speed)
function ke = compute_ke_from_meas(meas, sys, window_seconds)
if nargin < 3
    window_seconds = meas.t(end) - meas.t(1);
end
omega_ts = meas.omega_ts; % matrix [Nt x nGen]
% use delta w relative to mean pre-fault (assume first sample pre-fault near 0)
delta_omega = omega_ts - mean(omega_ts(1:min(3,size(omega_ts,1)),:),1);
% compute max delta omega per generator in window
domega_peak = max(abs(delta_omega),[],1)';
nGen = size(omega_ts,2);
ke = zeros(nGen,1);
for i=1:nGen
    M = sys.gen(i).M;
    ke(i) = 0.5 * M * (domega_peak(i))^2;
end
end
