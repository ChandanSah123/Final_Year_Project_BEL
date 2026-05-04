%% 1000Hz PSS/E PMU Data Verification Script
clc; clear; close all;

disp('--- Starting High-Resolution PMU Verification ---');
% 1. PYTHON BRIDGE SETUP
pe = pyenv('Version', 'C:\Users\Acer\AppData\Local\Programs\Python\Python312\python.exe', 'ExecutionMode', 'OutOfProcess');
if count(py.sys.path, pwd) == 0
    insert(py.sys.path, int32(0), pwd);
end

% Force reload the new engine to avoid cache issues
py.importlib.reload(py.importlib.import_module('lstm_predict'));

% 2. LOAD & PARSE NEW 1000Hz PMU DATA
loaded = load('data1.mat');
fields = fieldnames(loaded);
pmu_data = loaded.(fields{1}); 

% --- TIME PARAMETERS ---
dt_pmu = 0.001;           % NEW: PMU Resolution is now 1000 Hz
dt_ai  = 0.001;           % AI Native Resolution
t_clear = 1.16;            % Fault clears
t_target = 1.36;           % The exact time we want to verify (0.6s post-fault)

% 3. EXTRACT EXACTLY 200 ROWS FOR THE AI
disp('3. Extracting 200 rows (0.2s window) of 1000Hz data...');

% Calculate the exact row index for t = 1.2s
idx_start = round(t_clear / dt_pmu) + 1; 

% Slice exactly 200 rows (t=1.200 to t=1.399) 
% Columns: Angles (1-3), Speeds (4-6), Power (10-12)
ai_input_matrix = pmu_data(idx_start : idx_start + 199, [1:3, 4:6, 10:12]);

% Calculate exactly when this observation window ends so we can align the prediction
t_input_end = t_clear + (199 * dt_pmu); 

% Flatten and send to Python
flat_input = reshape(ai_input_matrix.', 1, []); 
py_input = py.list(flat_input);

% 4. RUN AI PREDICTION
disp('4. Predicting Post-Fault Trajectory...');
tic; 
py_output = py.lstm_predict.predict_next(py_input);
pred_flat = double(py_output);
time_taken = toc;

% 5. RECONSTRUCT & ALIGN TRAJECTORY
num_features = 6;
num_timesteps = length(pred_flat) / num_features;

if num_timesteps < 2
    disp('CRITICAL ERROR: Python returned a dummy array. Check Python code.');
    return;
end

% Rebuild the matrix [TimeSteps x Features]
pred_trajectory = reshape(pred_flat, num_features, num_timesteps).';

% Create the time vector for the AI's output (Column vector)
% It starts exactly one time-step after the input window ends
t_ai_output = linspace(t_input_end + dt_ai, t_input_end + (num_timesteps * dt_ai), num_timesteps).';

% Find the exact AI prediction at our target time (1.8s)
target_prediction = interp1(t_ai_output, pred_trajectory, t_target, 'linear');

pred_angles = target_prediction(1:3);
pred_speeds = target_prediction(4:6);

% 6. EXTRACT GROUND TRUTH FROM PSS/E
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