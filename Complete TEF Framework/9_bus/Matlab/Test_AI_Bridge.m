%% Minimal AI Deployment Test Script
clc; clear; close all;

disp('--- Starting Standalone AI Deployment Test ---');

% 1. Setup Python Bridge
% Force MATLAB to run Python in a completely isolated background process
pe = pyenv('Version', 'C:\Users\Acer\AppData\Local\Programs\Python\Python312\python.exe', 'ExecutionMode', 'OutOfProcess');

if count(py.sys.path, pwd) == 0
    insert(py.sys.path, int32(0), pwd);
end
disp('1. Python bridge connected (Out-Of-Process Mode Active).');

% 2. Generate Dummy Grid Data (200 steps, 9 features)
% In reality, this is: [Delta(1:3), Omega(1:3), Pe(1:3)]
disp('2. Generating 0.2 seconds of dummy grid data...');
raw_input_matrix = rand(200, 9); 

% 3. Format Data for Safe Transfer
% Flatten to 1D to completely avoid MATLAB-Python 2D matrix memory bugs
flat_input = reshape(raw_input_matrix.', 1, []); 
py_input = py.list(flat_input);
disp('3. Data flattened and formatted.');

% 4. Call Python and Predict
disp('4. Waking up Python AI Engine and predicting...');
tic; % Start a timer to see how fast the AI is
try
    % This calls the predict_next function inside lstm_predict.py
    py_output = py.lstm_predict1.predict_next(py_input);
    
    % 5. Receive and Parse Results
    pred_vals = double(py_output);
    time_taken = toc;
    
    disp(' ');
    disp('--- SUCCESS! AI PREDICTION RECEIVED ---');
    fprintf('Prediction Time: %.4f seconds\n', time_taken);
    disp('-----------------------------------------');
    fprintf('Predicted Generator Angles (Deg): [%.4f, %.4f, %.4f]\n', pred_vals(1), pred_vals(2), pred_vals(3));
    fprintf('Predicted Generator Speeds (PU):  [%.4f, %.4f, %.4f]\n', pred_vals(4), pred_vals(5), pred_vals(6));
    disp('-----------------------------------------');
    disp('If these numbers printed, your deployment is 100% working!');
    
catch ME
    disp(' ');
    disp('--- ERROR DURING PREDICTION ---');
    disp('The Python script failed to execute. Check the error message below:');
    disp(ME.message);
end