# lstm_predict.py
import os
import numpy as np
import joblib
import warnings
import traceback

# --- ANTI-FREEZE SETTINGS ---
os.environ['CUDA_VISIBLE_DEVICES'] = '-1' 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
warnings.filterwarnings('ignore')

from tensorflow.keras.models import load_model

# ==============================================================
# LOAD AI ASSETS ONCE
# ==============================================================
try:
    model = load_model('power_system_model.h5')
    scaler_X = joblib.load('scaler_X.pkl')
    scaler_y = joblib.load('scaler_y.pkl')
    print("Python AI Engine Loaded Successfully.")
except Exception as e:
    print(f"Error loading AI assets: {str(e)}")

def predict_next(flat_input_list):
    try:
        print("  -> Step 1 & 2: Parsing and Reshaping Input...")
        raw_data_flat = np.array(flat_input_list)
        raw_data_2d = raw_data_flat.reshape(200, 9)
        
        print(f"  -> Step 3: Scaling X... (Input Shape: {raw_data_2d.shape})")
        X_scaled = scaler_X.transform(raw_data_2d)
        
        print("  -> Step 4 & 5: Running AI Prediction...")
        X_lstm = X_scaled.reshape(1, 200, 9)
        pred_scaled = model.predict(X_lstm, verbose=0)
        
        print(f"  -> Step 5.5: Original AI Output Shape is {np.shape(pred_scaled)}")
        # ==========================================================
        # THE BRUTE-FORCE FIX
        # Crush any shape (3D, 4D, etc.) into a flat list, then grab  
        # the very last 6 numbers and force them into a 1x2D table.
        # ==========================================================
        pred_flat = np.array(pred_scaled).flatten()
        pred_2d = pred_flat[-6:].reshape(1, 6) 
        print(f"  -> Step 5.6: Forced to valid 2D Shape: {pred_2d.shape}")
        
        print("  -> Step 6: Inverse scaling to physical units...")
        final_pred = scaler_y.inverse_transform(pred_2d)
        
        print("  -> Step 7: Success! Returning to MATLAB.")
        return final_pred.flatten().tolist()
        
    except Exception as e:
        print("\n=== DETAILED PYTHON CRASH LOG ===")
        # This will print the exact line number of the crash
        traceback.print_exc()
        print("================================\n")
        return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


if __name__ == "__main__":
    print("Testing Python Script Standalone...")
    dummy_data = np.random.rand(1800).tolist()
    result = predict_next(dummy_data)
    print("\nPrediction Result:", result)