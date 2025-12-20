import psse3603
import psspy
import dyntools
import scipy.io
import numpy as np
import os

# ==============================================================================
# 1. SETUP
# ==============================================================================
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_NE\Python Script"
result_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_NE\Matlab Script"

# Update these filenames if they differ on your disk
raw_file = os.path.join(work_dir, "ieee39bus1.raw")
dyr_file = os.path.join(work_dir, "ieee39buscls.dyr")
out_file = os.path.join(result_dir, "NE.out")
mat_file = os.path.join(result_dir, "data1.mat")

# Initialize PSS/E
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
print("--- Starting Simulation (IEEE 39-Bus) ---")

# Simulation Parameters
fault_bus = 22
fault_time = 1.0
clear_time = 1.3  # Note: 0.3s fault is very long, might be unstable
end_time = 4.0

# Load Case
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

# Solver Settings
psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.01, _f, _f, _f, _f, _f])
psspy.fnsl([0,0,0,1,1,0,99,0])

# --- CHANNEL SETUP ---
# We use wildcard [-1, -1, -1, 1, TYPE, 0] to select ALL 10 generators automatically.
# The order here defines the column order in data1.mat

# 1. ANGLE (Code 1) -> Cols 1-10
psspy.chsb(0, 1, [-1, -1, -1, 1, 1, 0])

# 2. MECHANICAL POWER (Code 3) -> Cols 11-20
# (Your previous code used 3, which is correct for Pmech)
psspy.chsb(0, 1, [-1, -1, -1, 1, 6, 0])

# 3. ELECTRICAL POWER (Code 2) -> Cols 21-30
psspy.chsb(0, 1, [-1, -1, -1, 1, 2, 0])

# 4. SPEED (Code 7) -> Cols 31-40
psspy.chsb(0, 1, [-1, -1, -1, 1, 7, 0])

# --- CONVERT & RUN ---
print("Converting Loads...")
psspy.cong(0)
psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])

print("Running Simulation...")
psspy.strt(0, out_file)
psspy.run(0, fault_time, 0, 1, 0)

print(f"Applying Fault at Bus {fault_bus}...")
psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])

psspy.run(0, clear_time, 0, 1, 0)

print("Clearing Fault...")
psspy.dist_clear_fault(1)

psspy.run(0, end_time, 0, 1, 0)
print("Simulation Complete.")

# ==============================================================================
# 3. DATA EXTRACTION (Strictly data1.mat)
# ==============================================================================
print("Extracting Data...")
chnf_obj = dyntools.CHNF(out_file)
short_title, chanid, chandata = chnf_obj.get_data()
t_data = chandata['time']

# Filter only integer keys (channel IDs) and sort them (1 to 40)
sorted_keys = sorted([k for k in chanid.keys() if isinstance(k, int)])

data_columns = []

# Loop to build the matrix
for key in sorted_keys:
    y_data = chandata[key]
    min_len = min(len(t_data), len(y_data))
    data_columns.append(y_data[:min_len])

# Stack columns horizontally
if data_columns:
    final_matrix = np.column_stack(data_columns)
    
    # Save to .mat
    scipy.io.savemat(mat_file, {"data": final_matrix})
    
    print(f"SUCCESS: {mat_file} saved.")
    print(f"Matrix Shape: {final_matrix.shape}")
    print("Format:")
    print("  Cols 1-10:   Rotor Angle")
    print("  Cols 11-20:  Mechanical Power")
    print("  Cols 21-30:  Electrical Power")
    print("  Cols 31-40:  Rotor Speed")
else:
    print("Error: No data extracted.")