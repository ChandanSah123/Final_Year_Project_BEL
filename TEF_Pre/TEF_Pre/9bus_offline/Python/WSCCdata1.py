try:
    import psse3603  # type: ignore
    import psspy     # type: ignore
    import dyntools  # type: ignore
    pss_available = True
except Exception:
    pss_available = False
import scipy.io
import numpy as np
import os
import config
import pandas as pd
# ==============================================================================
# 1. SETUP
# ==============================================================================
work_dir = config.WORK_DIR
result_dir = config.RESULT_DIR
raw_file = config.RAW_FILE
dyr_file = config.DYR_FILE
out_file = config.OUT_FILE.replace('IEEE9_Combined.out','IEEE9.out')
mat_file = config.MAT_FILE_DATA
Tcl_file=config.TIME_FILE

# Initialize PSS/E (safe)
try:
    _i = psspy.getdefaultint()
    _f = psspy.getdefaultreal()
    psspy.psseinit(50)
except Exception:
    _i = 0
    _f = 0.0

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
print("--- Starting Simulation for IEEE 9-Bus ---")
print("--- Starting Simulation for IEEE 39-Bus ---")

# --- Simulation Settings ---
fault_bus = config.FAULT_BUS
from_bus = config.TRIP_LINE_FROM
to_bus = config.TRIP_LINE_TO
line_id = config.LINE_ID

fault_time = config.T_FAULT_START
t_cl = config.T_CLEAR
clear_time = fault_time + t_cl
end_time = config.T_END

# --- Load Case ---
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

# --- Solve Power Flow ---
psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
psspy.fnsl([0,0,0,1,1,0,99,0])

# ==============================================================================
# 3. CHANNEL SETUP (10 GENERATORS)
# ==============================================================================
psspy.delete_all_plot_channels()

# 1. ANGLE (Code 1) -> Becomes Channels 1-10
psspy.chsb(0, 1, [-1, -1, -1, 1, 1, 0])

# 2. SPEED (Code 7) -> Becomes Channels 11-20
psspy.chsb(0, 1, [-1, -1, -1, 1, 7, 0])

# 3. MECHANICAL POWER (Code 6) -> Becomes Channels 21-30
psspy.chsb(0, 1, [-1, -1, -1, 1, 6, 0])

# 4. ELECTRICAL POWER (Code 2) -> Becomes Channels 31-40
psspy.chsb(0, 1, [-1, -1, -1, 1, 2, 0])

# 5. REACTIVE POWER (Code 4) -> Becomes Channels 41-50
psspy.chsb(0, 1, [-1, -1, -1, 1, 4, 0])

# ==============================================================================
# 4. CONVERT & RUN
# ==============================================================================
print("Converting Network...")
psspy.cong(0) 
psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.fact()
psspy.tysl(0)

# Run
ierr = psspy.strt(0, out_file)
if ierr != 0:
    print(f"Error starting simulation: {ierr}")
    exit()

psspy.run(0, fault_time, 0, 1, 0)

print(f"Applying Fault at Bus {fault_bus}...")
psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])

psspy.run(0, clear_time, 0, 1, 0)

print(f"Tripping Line {from_bus}-{to_bus}...")
#psspy.dist_branch_trip(from_bus, to_bus, line_id)
psspy.dist_clear_fault(1)

psspy.run(0, end_time, 0, 1, 0)

# ==============================================================================
# 5. EXTRACTION & EXPORT
# ==============================================================================
print("Extracting Data...")
try:
    chnf_obj = dyntools.CHNF(out_file)
    short_title, chanid, chandata = chnf_obj.get_data()
except Exception as e:
    print(f"Error extracting data: {e}")
    exit()

t_data = chandata['time']

# Filter only integer keys (channel IDs 1-50) and sort them
sorted_keys = sorted([k for k in chanid.keys() if isinstance(k, int)])

data_columns = []

# 1. Collect standard data channels (Cols 1-50)
for key in sorted_keys:
    y_data = chandata[key]
    # Ensure length matches time (truncate if needed)
    min_len = min(len(t_data), len(y_data))
    data_columns.append(y_data[:min_len])

if data_columns:
    # 2. Append Time as the LAST column (Column 51)
    # We use min_len from the loop to ensure time vector matches data length
    # (Assuming all channels have same length, which they do in PSS/E)
    data_len = len(data_columns[0])
    time_col = t_data[:data_len]
    data_columns.append(time_col)

    # Stack columns horizontally
    final_matrix = np.column_stack(data_columns)
    
    # Save to .mat
    scipy.io.savemat(mat_file, {"data": final_matrix})

    # --- VERIFY MAPPING ---
    print("\n--- VERIFYING CHANNEL MAPPING ---")
    for i in range(1, 4):
        if i in chanid:
            print(f"Column {i}: {chanid[i]}")

    print("-" * 30)
    print(f"SUCCESS: {mat_file} saved.")
    print("DATA STRUCTURE MAPPING:")
    print(f"Matrix Shape: {final_matrix.shape}")
    print("Columns 1-3:   Rotor Angle (delta)")
    print("Columns 4-6:   Speed (omega)")
    print("Columns 7-9:   Mechanical Power (Pm)")
    print("Columns 10-12: Electrical Power (Pe)")
    print("Columns 13-15: Reactive Power (Q)")
    print("Column  16:    Simulation Time (t)")
else:
    print("Error: No data extracted.")

    # ... (existing code for saving data matrix) ...

# ==============================================================================
# 6. EXPORT CLEARING TIME
# ==============================================================================
try:
    # Create a dictionary with the variable name you want in MATLAB
    # We wrap t_cl in a list [t_cl] or float(t_cl) to ensure compatible format
    tcl_data = {"t_cl": float(t_cl)}
    
    # Save to the file path defined in config.TIME_FILE
    scipy.io.savemat(Tcl_file, tcl_data)
    
    print(f"SUCCESS: Clearing time ({t_cl} s) saved to {Tcl_file}")

except Exception as e:
    print(f"Error saving clearing time: {e}")