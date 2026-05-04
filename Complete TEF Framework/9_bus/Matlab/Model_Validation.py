import os
import numpy as np
import joblib
import scipy.io as sio
import traceback
from tensorflow.keras.models import load_model

# --- Environment Setup ---
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

def run_full_validation():
    # 1. SETTINGS & PATHS
    mat_file = 'data1.mat'
    model_file = 'power_system_model.h5'
    scaler_x_file = 'scaler_X.pkl'
    scaler_y_file = 'scaler_y.pkl'
    
    t_clear_idx = 1200   # Start of window (1.2s at 1000Hz)
    window_size = 200    # 0.2s input window
    target_idx = 1800    # Validation point (1.8s at 1000Hz)
    
    # Column indices for X (Input): Angles(0-2), Speeds(3-5), Power(9-11)
    cols_x = [0, 1, 2, 3, 4, 5, 9, 10, 11]
    
    try:
        # 2. LOAD DATA & ASSETS
        print(f"--- Loading {mat_file} ---")
        mat_contents = sio.loadmat(mat_file)
        data_key = [k for k in mat_contents.keys() if not k.startswith('_')][0]
        pmu_data = mat_contents[data_key]
        
        print("--- Loading AI Assets ---")
        model = load_model(model_file, compile=False)
        scaler_X = joblib.load(scaler_x_file)
        scaler_y = joblib.load(scaler_y_file)
        
        # 3. PREPARE INPUT WINDOW
        input_raw = pmu_data[t_clear_idx : t_clear_idx + window_size, cols_x]
        X_scaled = scaler_X.transform(input_raw)
        X_lstm = X_scaled.reshape(1, window_size, 9)
        
        # 4. RUN PREDICTION
        print("--- Running Inference ---")
        pred_scaled = model.predict(X_lstm, verbose=0)
        
        # 5. INVERSE SCALE & EXTRACT TARGET
        # Output is (1, 800, 6). Index 600 in output corresponds to T=1.8s 
        # (since output starts at 1.2s + 0.2s window = 1.4s)
        final_trajectory = scaler_y.inverse_transform(pred_scaled[0])
        pred_at_target = final_trajectory[400] # Adjusting: 1.4s to 1.8s is 400 steps
        
        # 6. GET GROUND TRUTH
        actual_at_target = pmu_data[target_idx, 0:6]
        
        # 7. FORMATTED RESULTS
        print("\n" + "="*60)
        print(f"   VALIDATION RESULTS AT T = 1.8s (Row {target_idx})")
        print("="*60)
        
        print(f"{'Metric':<20} | {'Actual (PSS/E)':<15} | {'Predicted (AI)':<15} | {'Abs Error'}")
        print("-" * 60)
        
        labels = ['Angle 1 (Deg)', 'Angle 2 (Deg)', 'Angle 3 (Deg)', 
                  'Speed 1 (PU)', 'Speed 2 (PU)', 'Speed 3 (PU)']
        
        for i in range(6):
            actual = actual_at_target[i]
            predicted = pred_at_target[i]
            error = abs(actual - predicted)
            print(f"{labels[i]:<20} | {actual:<15.4f} | {predicted:<15.4f} | {error:.4f}")
            
        print("="*60)
        print(f"Trajectory sequence length: {len(final_trajectory)} steps")

    except Exception as e:
        print("\n!!! ERROR OCCURRED !!!")
        traceback.print_exc()

if __name__ == "__main__":
    run_full_validation()