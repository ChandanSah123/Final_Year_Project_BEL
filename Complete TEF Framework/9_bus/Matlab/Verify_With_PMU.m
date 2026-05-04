%% Real PMU Data Verification Script (Time-Synchronized)
clc; clear; close all;

disp('--- Starting Synchronized PMU Verification ---');

% 1. PYTHON BRIDGE SETUP
pe = pyenv('Version', 'C:\Users\Acer\AppData\Local\Programs\Python\Python312\python.exe', 'ExecutionMode', 'OutOfProcess');
if count(py.sys.path, pwd) == 0
    insert(py.sys.path, int32(0), pwd);
end

% ==============================================================
% THE FIX: Force MATLAB to clear its cache and reload the script
% ==============================================================
py.importlib.reload(py.importlib.import_module('lstm_predict'));
disp('Python module forcefully reloaded!');

% 2. LOAD & PARSE PMU DATA
loaded = load('data1.mat');
fields = fieldnames(loaded);
pmu_data = loaded.(fields{1}); 

% --- TIME PARAMETERS ---
dt_pmu = 0.01;            % PMU Resolution (100 Hz)
dt_ai  = 0.001;           % AI Native Resolution (1000 Hz)
t_clear = 1.2;            % Fault clears
t_input_end = 1.4;        % End of observation window
t_target = 1.8;           % The exact time we want to verify (0.6s post-fault)

% Calculate PMU indices for the 0.2s input window
idx_start = round(t_clear / dt_pmu) + 1; 
idx_end   = round(t_input_end / dt_pmu) + 1; 

% 3. EXTRACT AND RESAMPLE FEATURES
disp('3. Resampling PMU Data (100Hz -> 1000Hz)...');
% Grab the exact 0.2 second window (approx 21 rows)
raw_pmu_window = pmu_data(idx_start:idx_end, [1:3, 4:6, 10:12]);

% Create time vectors for interpolation
t_pmu_vector = linspace(t_clear, t_input_end, size(raw_pmu_window, 1));
t_ai_vector  = linspace(t_clear, t_input_end, 200);

% Up-sample to 200 points to match AI training
ai_input_matrix = interp1(t_pmu_vector, raw_pmu_window, t_ai_vector, 'spline');

% Flatten and send to Python
flat_input = reshape(ai_input_matrix.', 1, []); 
py_input = py.list(flat_input);

% 4. RUN AI PREDICTION
disp('4. Predicting Post-Fault Trajectory...');
tic; 
py_output = py.lstm_predict.predict_next(py_input);
pred_flat = double(py_output);
time_taken = toc;

% ==============================================================
% 5. RECONSTRUCT & ALIGN TRAJECTORY
% ==============================================================
num_features = 6;
num_timesteps = length(pred_flat) / num_features;

% --- EXPOSE THE DATA SIZE ---
fprintf('\n--- DIAGNOSTICS ---\n');
fprintf('Data points received from Python: %d (Should be 4800)\n', length(pred_flat));
fprintf('Calculated AI timesteps: %d (Should be 800)\n', num_timesteps);
disp('-------------------');

if num_timesteps < 2
    disp('CRITICAL ERROR: Python returned the 6-number dummy array!');
    disp('This means the Python script crashed internally. Check your Python code.');
    return; % Stop the script safely
end

% Rebuild the matrix [TimeSteps x Features]
pred_trajectory = reshape(pred_flat, num_features, num_timesteps).';

% Create the time vector as a COLUMN vector (.') to perfectly match the matrix rows
t_ai_output = linspace(t_input_end, t_input_end + (num_timesteps * dt_ai), num_timesteps).';

% Find the exact AI prediction at our target time (1.8s)
target_prediction = interp1(t_ai_output, pred_trajectory, t_target, 'linear');

pred_angles = target_prediction(1:3);
pred_speeds = target_prediction(4:6);

% 6. EXTRACT GROUND TRUTH
idx_target_pmu = round(t_target / dt_pmu) + 1;
actual_angles = pmu_data(idx_target_pmu, 1:3);
actual_speeds = pmu_data(idx_target_pmu, 4:6);

error_angles = abs(actual_angles - pred_angles);
error_speeds = abs(actual_speeds - pred_speeds);

% 7. DISPLAY RESULTS
disp(' ');
disp('====================================================');
fprintf('    AI VERIFICATION RESULTS (Target T = %.2fs)      \n', t_target);
disp('====================================================');
fprintf('Inference Time: %.4f seconds\n\n', time_taken);

disp('--- GENERATOR ANGLES (Deg) ---');
fprintf('Predicted: [%.4f, %.4f, %.4f]\n', pred_angles(1), pred_angles(2), pred_angles(3));
fprintf('Actual:    [%.4f, %.4f, %.4f]\n', actual_angles(1), actual_angles(2), actual_angles(3));
fprintf('Error:     [%.4f, %.4f, %.4f]\n\n', error_angles(1), error_angles(2), error_angles(3));

disp('--- GENERATOR SPEEDS (PU) ---');
fprintf('Predicted: [%.4f, %.4f, %.4f]\n', pred_speeds(1), pred_speeds(2), pred_speeds(3));
fprintf('Actual:    [%.4f, %.4f, %.4f]\n', actual_speeds(1), actual_speeds(2), actual_speeds(3));
fprintf('Error:     [%.4f, %.4f, %.4f]\n', error_speeds(1), error_speeds(2), error_speeds(3));
disp('====================================================');

%%Predicted: [109.7839, 132.7101, 143.2086]
%Actual:    [57.5975, 421.0154, 471.6301]
%Error:     [52.1863, 288.3053, 328.4214]

%--- GENERATOR SPEEDS (PU) ---
%Predicted: [0.0077, 0.0060, 0.0106]
%Actual:    [0.0001, 0.0578, 0.0489]
%Error:     [0.0076, 0.0518, 0.0383]
%====================================================