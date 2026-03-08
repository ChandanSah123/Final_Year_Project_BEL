function [MOD_sort_data, t_pebs_idx, PEBS_Crossing] = detect_pebs_crossing(fth, t_sim, fault_start_time, npts, th, w_tilde)
% DETECT_PEBS_CROSSING (Modified for Grazing/Near-Miss Detection)
% Finds PEBS crossing (Standard Zero-Cross OR Grazing/Local Peak).
%
% Output:
%   MOD_sort_data : Sorted by Angle [ID, Angle, Speed]
%   t_pebs_idx    : Index of crossing
%   PEBS_Crossing : Unsorted [ID, Angle, Speed] (1 to g)

    % Initialize outputs
    MOD_sort_data = [];
    PEBS_Crossing = [];
    t_pebs_idx = 0;
    
    % PARAMETER: Sensitivity Threshold
    % If fth gets closer than this to 0 (e.g., -0.05) and turns back down, 
    % we count it as a crossing (Potential Energy Ridge reached).
    grazing_threshold = -0.05; 

    % 1. DEFINE START POINT
    idx_start = find(t_sim >= fault_start_time, 1);
    if isempty(idx_start)
        idx_start = 2; 
    else
        idx_start = idx_start + 1; 
    end
    
    t = idx_start + 1;
    found_crossing = false;
    detection_type = ''; % To store if it was 'crossing' or 'grazing'

    % 2. Find the Crossing
    % Note: We loop until npts-1 so we can safely check t+1 for the peak detection
    while t < npts 
        curr_val = fth(t);
        prev_val = fth(t-1);
        next_val = fth(t+1);

        % --- CHECK A: Standard Zero Crossing ---
        % The curve goes from negative to positive
        if (curr_val > 0 && prev_val <= 0)
            t_pebs_idx = t;
            found_crossing = true;
            detection_type = 'Standard Zero Crossing';
            break;
        end

        % --- CHECK B: Grazing / Local Peak ---
        % The curve stays negative but gets very close to zero and turns down.
        % 1. Still negative (curr_val < 0)
        % 2. Higher than threshold (curr_val > -0.05)
        % 3. Is a local peak (higher than previous and next points)
        if (curr_val < 0) && (curr_val > grazing_threshold) && ...
           (curr_val >= prev_val) && (curr_val > next_val)
       
            t_pebs_idx = t;
            found_crossing = true;
            detection_type = 'Grazing (Local Peak)';
            break;
        end
        
        t = t + 1;
    end

    if found_crossing
        % 3. Interpolate State at Crossing
        
        if strcmp(detection_type, 'Standard Zero Crossing')
            % Standard Logic: Interpolate between t-1 and t
            y1 = fth(t_pebs_idx-1);
            y2 = fth(t_pebs_idx);
            fraction = abs(y1) / (abs(y1) + abs(y2));
            
            % Interpolate Time
            crossing_time = t_sim(t_pebs_idx-1) + fraction * (t_sim(t_pebs_idx) - t_sim(t_pebs_idx-1));
            
            % Interpolate Angles
            raw_angles_prev = th(t_pebs_idx-1, :);
            raw_angles_curr = th(t_pebs_idx, :);
            raw_angles = raw_angles_prev + fraction * (raw_angles_curr - raw_angles_prev);
            
            % Interpolate Speeds
            raw_speeds_prev = w_tilde(t_pebs_idx-1, :);
            raw_speeds_curr = w_tilde(t_pebs_idx, :);
            raw_speeds = raw_speeds_prev + fraction * (raw_speeds_curr - raw_speeds_prev);
            
        else
            % Grazing Logic: We are AT the peak. No interpolation between t-1 and t needed.
            % We just take the values at the peak index (t).
            crossing_time = t_sim(t_pebs_idx);
            raw_angles = th(t_pebs_idx, :);
            raw_speeds = w_tilde(t_pebs_idx, :);
        end
        
        gen_ids = 1:length(raw_angles);   
        
        % 4. Create Unsorted Matrix (PEBS_Crossing)
        PEBS_Crossing = [gen_ids', raw_angles', raw_speeds'];
        
        % 5. Sort by Angle for MOD Analysis
        % Note: Usually sorting by Angle magnitude or Speed (depending on method)
        % Your original code sorts by Angle value descending.
        [sorted_angles, sort_idx] = sort(raw_angles, 'descend');
        sorted_speeds = raw_speeds(sort_idx);
        sorted_ids    = gen_ids(sort_idx);
        
        MOD_sort_data = [sorted_ids', sorted_angles', sorted_speeds'];
        
        % 6. Save and Display
        save('C_MOD.mat', 'MOD_sort_data', 'PEBS_Crossing');
        
        fprintf('Success! PEBS Crossing found via [%s] at t = %.4f s (Index %d).\n', detection_type, crossing_time, t_pebs_idx);
        
        disp('Variable "PEBS_Crossing" (Unsorted) content:');
        disp(PEBS_Crossing);
        
        disp('Variable "MOD_sort_data" (Sorted) content:');
        disp(MOD_sort_data);
        
    else
        disp('No PEBS crossing found (Standard or Grazing).');
    end
end