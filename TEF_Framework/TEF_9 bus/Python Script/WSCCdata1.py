import psse3603
import psspy
import dyntools
import scipy.io
import numpy as np
import os

# ==============================================================================
# 1. SETUP
# ==============================================================================
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_9 bus\Python Script"
result_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_9 bus\Matlab Script"
raw_file = os.path.join(work_dir, "IEEE9bus.raw")
dyr_file = os.path.join(work_dir, "ieee9bus.dyr")
out_file = os.path.join(result_dir, "IEEE9.out")
mat_file = os.path.join(result_dir, "data1.mat") # Saving directly as data1.mat

# Initialize PSS/E
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
print("--- Starting Simulation for IEEE 9-Bus ---")

# Simulation Settings
fault_bus = 7
from_bus = 7
to_bus = 5
line_id = '1'
fault_time = 1.0
clear_time = 1.2
end_time = 4.0

# Load Case
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

# Solve Power Flow
psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.01, _f, _f, _f, _f, _f])
psspy.fnsl([0,0,0,1,1,0,99,0])

# --- CHANNEL SETUP (CRITICAL STEP) ---
# The order here determines the column order in the MAT file.
# We use -1, -1, -1 to select "All Generators" automatically.

# 1. ANGLE (Code 1) -> Becomes Channels 1, 2, 3
psspy.chsb(0, 1, [-1, -1, -1, 1, 1, 0])

# 2. MECHANICAL POWER (Code 3) -> Becomes Channels 4, 5, 6
# (Note: Previous code used 6 for Voltage. We changed it to 3 for Pmech)
psspy.chsb(0, 1, [-1, -1, -1, 1, 6, 0])

# 3. ELECTRICAL POWER (Code 2) -> Becomes Channels 7, 8, 9
psspy.chsb(0, 1, [-1, -1, -1, 1, 2, 0])

# 4. SPEED (Code 7) -> Becomes Channels 10, 11, 12
psspy.chsb(0, 1, [-1, -1, -1, 1, 7, 0])

# --- CONVERT & RUN ---
psspy.cong(0) 
psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])

psspy.strt(0, out_file)
psspy.run(0, fault_time, 0, 1, 0)

print(f"Fault at Bus {fault_bus}...")
psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])

psspy.run(0, clear_time, 0, 1, 0)

print(f"Tripping Line {from_bus}-{to_bus}...")
psspy.dist_branch_trip(from_bus, to_bus, line_id)
psspy.dist_clear_fault(1)

psspy.run(0, end_time, 0, 1, 0)

# ==============================================================================
# 3. EXTRACTION & EXPORT (DATA1.MAT ONLY)
# ==============================================================================
print("Extracting Data...")
chnf_obj = dyntools.CHNF(out_file)
short_title, chanid, chandata = chnf_obj.get_data()
t_data = chandata['time']

# Filter only integer keys (channel IDs) and sort them (1 to 12)
sorted_keys = sorted([k for k in chanid.keys() if isinstance(k, int)])

data_columns = []

# Loop 1-12 to build the matrix
for key in sorted_keys:
    y_data = chandata[key]
    # Ensure all columns are same length as time
    min_len = min(len(t_data), len(y_data))
    data_columns.append(y_data[:min_len])

# Stack columns horizontally
# Shape will be [TimeSteps x 12]
if data_columns:
    final_matrix = np.column_stack(data_columns)
    
    # Save to .mat
    scipy.io.savemat(mat_file, {"data": final_matrix})
    
    print(f"SUCCESS: {mat_file} saved.")
    print(f"Matrix Shape: {final_matrix.shape}")
    print("Columns 1-3: Angle")
    print("Columns 4-6: P_mech")
    print("Columns 7-9: P_elec")
    print("Columns 10-12: Speed")
else:
    print("Error: No data extracted.")