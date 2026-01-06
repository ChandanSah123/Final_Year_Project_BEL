import os
import sys
import textwrap
import scipy.io
import numpy as np
import config

# PSS/E Initialization
try:
    import psse3603 # type: ignore
    import psspy    # type: ignore
    import dyntools # type: ignore
except ImportError:
    sys.path.append(r"C:\Program Files\PTI\PSSE35\PSSPY37")
    sys.path.append(r"C:\Program Files\PTI\PSSE35\PSSBIN")
    import psspy
    import dyntools

_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(1000)

# ==============================================================================
# 1. CONFIGURATION: BUS FAULT (NO TRIP)
# ==============================================================================
raw_file = config.RAW_FILE
dyr_file = config.DYR_FILE
# Save as data1.mat (or modify name if you want distinct files)
mat_file = config.MAT_FILE_DATA 
out_file = config.OUT_FILE.replace('.out', '_BusFault.out')

fault_bus = config.FAULT_BUS
t_fault   = config.T_FAULT_START
t_clear   = config.T_FAULT_START + config.T_CLEAR
t_end     = config.T_END

print("--- Starting Simulation: Bus Fault (No Trip) ---")

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
# Load Case
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

# Solve Power Flow
psspy.fnsl([0,0,0,1,1,0,99,0])

# Define Channels (Robust Method)
psspy.delete_all_plot_channels()
# 1. ANGLE (Code 1) -> Cols 1-3
psspy.chsb(0, 1, [-1, -1, -1, 1, 1, 0])
# 2. SPEED (Code 7) -> Cols 4-6
psspy.chsb(0, 1, [-1, -1, -1, 1, 7, 0])
# 3. MECHANICAL POWER (Code 6) -> Cols 7-9
psspy.chsb(0, 1, [-1, -1, -1, 1, 6, 0])
# 4. ELECTRICAL POWER (Code 2) -> Cols 10-12
psspy.chsb(0, 1, [-1, -1, -1, 1, 2, 0])
# 5. REACTIVE POWER (Code 4) -> Cols 13-15
psspy.chsb(0, 1, [-1, -1, -1, 1, 4, 0])

# Convert Network
psspy.cong(0)
psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.fact()
psspy.tysl(0)

# Run Simulation
psspy.strt(0, out_file)
psspy.run(0, t_fault, 0, 1, 0)

print(f"Applying Fault at Bus {fault_bus}...")
psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])

psspy.run(0, t_clear, 0, 1, 0)

print("Clearing Fault (No Trip)...")
psspy.dist_clear_fault(1)

psspy.run(0, t_end, 0, 1, 0)

# ==============================================================================
# 3. EXTRACTION & EXPORT
# ==============================================================================
print("Extracting Data...")
try:
    chnf_obj = dyntools.CHNF(out_file)
    short_title, chanid, chandata = chnf_obj.get_data()
except Exception as e:
    print(f"Error extracting data: {e}")
    exit()

t_data = chandata['time']
# Filter only integer keys (channel IDs) and sort them
sorted_keys = sorted([k for k in chanid.keys() if isinstance(k, int)])
data_columns = []

for key in sorted_keys:
    y_data = chandata[key]
    min_len = min(len(t_data), len(y_data))
    data_columns.append(y_data[:min_len])

if data_columns:
    # Append Time as Last Column
    data_len = len(data_columns[0])
    time_col = t_data[:data_len]
    data_columns.append(time_col)

    final_matrix = np.column_stack(data_columns)
    
    scipy.io.savemat(mat_file, {"data": final_matrix})
    
    print("-" * 30)
    print(f"SUCCESS: {mat_file} saved.")
    print(f"Matrix Shape: {final_matrix.shape}")
    print("Columns 1-3:   Rotor Angle")
    print("Columns 4-6:   Speed")
    print("Columns 7-9:   Pm")
    print("Columns 10-12: Pe")
    print("Columns 13-15: Q")
    print("Column  16:    Time")
else:
    print("Error: No data extracted.")