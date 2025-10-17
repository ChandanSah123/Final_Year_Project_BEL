% TEF_offline_build.m
% Builds the offline DB: for each OP and fault, simulate fault, compute KE fingerprints,
% compute CUEP/Vcr (BCU heavy), store in DB struct.
function DB = TEF_offline_build(sys, opts)
% opts.operating_points: cell array of OP names or parameter structs
% opts.fault_list: cell array of fault identifiers
% opts.offline_dt: float seconds - window to compute KE fingerprint
DB = struct();
entry_idx = 0;
fprintf('Starting offline DB build...\n');
for op_i = 1:numel(opts.operating_points)
    OP = opts.operating_points{op_i};
    % TODO: if multiple OPs, vary loads/generators accordingly (power flow)
    for f = 1:numel(opts.fault_list)
        faultLoc = opts.fault_list{f};
        % create a standard faultInfo for offline (choose representative durations)
        faultInfo.location = faultLoc;
        faultInfo.t_fault = 0.0;
        faultInfo.t_clear = 0.08; % example
        simOpts.tfinal = 1.0;
        simOpts.dt = 1e-3;
        [traj, meas] = simulate_fault(sys, faultInfo, simOpts);
        % compute KE fingerprint over first opts.offline_dt seconds
        ke = compute_ke_from_meas(meas, sys, opts.offline_dt);
        KE_norm = ke / sum(ke + eps);
        % determine nominal MOD from full simulation (who diverges in long run)
        MOD = detect_MOD_from_traj(traj, sys);
        % compute post-fault SEP
        Ypost = util_Ybus_update(sys.Ybus, faultInfo);
        theta_s = compute_postfault_SEP(Ypost, sys);
        % compute CUEP and Vcr (offline heavy)
        try
            [theta_u, Vcr] = compute_CUEP_BCU(sys, Ypost, theta_s, MOD);
        catch
            warning('BCU failed; using corner-point initial guess and refine');
            theta_u = compute_CUEP_via_corner_and_refine(MOD, theta_s, sys, Ypost);
            Vcr = compute_energy(theta_u, zeros(size(theta_u)), theta_s, sys);
        end
        % store
        entry_idx = entry_idx + 1;
        DB(entry_idx).fault = faultLoc;
        DB(entry_idx).op = OP;
        DB(entry_idx).KE_norm = KE_norm;
        DB(entry_idx).MOD = MOD;
        DB(entry_idx).theta_u = theta_u;
        DB(entry_idx).Vcr = Vcr;
        DB(entry_idx).theta_s = theta_s;
        DB(entry_idx).topologyHash = topology_hash(Ypost);
        fprintf('DB entry %d: fault=%s, OP=%s, MOD=%s\n', entry_idx, faultLoc, OP, mat2str(MOD));
    end
end
fprintf('Offline DB build complete. Entries: %d\n', numel(DB));
end

% Helper: very simple topology hashing (replaceable)
function h = topology_hash(Y)
h = sprintf('%08x', sum(abs(Y(:))) * 1e6);
end
