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
        raw_data_flat = np.array(flat_input_list)
        raw_data_2d = raw_data_flat.reshape(200, 9)
        
        X_scaled = scaler_X.transform(raw_data_2d)
        X_lstm = X_scaled.reshape(1, 200, 9)
        
        pred_scaled = model.predict(X_lstm, verbose=0)
        
        # ==========================================================
        # NEW FIX: Return the ENTIRE trajectory (800 x 6)
        # ==========================================================
        if len(pred_scaled.shape) == 3:
            pred_2d = pred_scaled[0] # Grab the full 800x6 matrix
        else:
            pred_2d = pred_scaled
            
        # Inverse transform the entire trajectory back to physical units
        final_pred = scaler_y.inverse_transform(pred_2d)
        
        # Flatten the massive matrix and send it all to MATLAB
        return final_pred.flatten().tolist()
        
    except Exception as e:
        print("\n=== DETAILED PYTHON CRASH LOG ===")
        traceback.print_exc()
        print("================================\n")
        # Return a safe dummy array if something crashes
        return [0.0] * 6


if __name__ == "__main__":
    print("Testing Python Script Standalone...")
    dummy_data = np.random.rand(1800).tolist()
    result = predict_next(dummy_data)
    print(f"\nPrediction Result Length: {len(result)} numbers (Should be 4800)")