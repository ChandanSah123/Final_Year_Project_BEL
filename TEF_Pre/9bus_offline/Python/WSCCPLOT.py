try:
    import psse3603  # type: ignore
    import psspy     # type: ignore
    import dyntools  # type: ignore
    pss_available = True
except Exception:
    pss_available = False
import matplotlib.pyplot as plt
import pandas as pd
import scipy.io
import numpy as np
import os
import re
import config

# ==============================================================================
# 1. SETUP & CONFIGURATION
# ==============================================================================
work_dir = config.WORK_DIR
result_dir = config.RESULT_DIR
raw_file = config.RAW_FILE
dyr_file = config.DYR_FILE
out_file = config.OUT_FILE.replace('IEEE9_Combined.out','IEEE9.out')

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

# Simulation Parameters

fault_bus = config.FAULT_BUS
from_bus = config.TRIP_LINE_FROM
to_bus = config.TRIP_LINE_TO
line_id = config.LINE_ID
fault_time = config.T_FAULT_START
t_cl = config.T_CLEAR
clear_time = fault_time + t_cl
end_time = config.T_END

# Load Case
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

# Solver Settings
psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
psspy.fnsl([0,0,0,1,1,0,99,0])

# --- DEFINE CHANNELS (UPDATED) ---
# NOTE: Assumes Generators are at Bus 1, 2, and 3 with ID '1'

#Channels

psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,6,0])
psspy.chsb(0,1,[-1,-1,-1,1,2,0])
psspy.chsb(0,1,[-1,-1,-1,1,7,0])

# --- CONVERT NETWORK ---
print("Converting Network...")
psspy.cong(0) 

# Load Conversion (3 Steps)
psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])

# Run Simulation
print("Initializing State (STRT)...")
ierr_strt = psspy.strt(0, out_file)

if ierr_strt > 0:
    print(f"ERROR: Simulation failed to start. Error code: {ierr_strt}")
else:
    psspy.run(0, fault_time, 0, 1, 0)

    print(f"Applying Fault at Bus {fault_bus}...")
    psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])

    psspy.run(0, clear_time, 0, 1, 0)
    print("Clearing Fault by TRIPPING LINE 7-5...")
    psspy.dist_branch_trip(from_bus, to_bus, line_id)   # Use correct line ID

    print("Clearing Fault...")
    psspy.dist_clear_fault(1)

    psspy.run(0, end_time, 0, 1, 0)
    print("Simulation Complete.")

# ==============================================================================
# 3. PLOTTING
# ==============================================================================

if not os.path.exists(out_file):
    print(f"CRITICAL ERROR: Output file not found at {out_file}")
else:
    chnf_obj = dyntools.CHNF(out_file)
    short_title, chanid, chandata = chnf_obj.get_data()
    time_data = chandata['time']

    print("Generating Plots...")

    # Plot 1: Rotor Angles (Channels 1, 2, 3)
    plt.figure(figsize=(10, 6))
    for ch in range(1, 4):
        if ch in chandata:
            clean_label = str(chanid[ch]).strip()
            plt.plot(time_data, chandata[ch], label=f"Ch{ch}: {clean_label}")
    plt.title(f"Rotor Angles (CT = {t_cl:.4f}s)")
    plt.ylabel("Angle (Deg)")
    plt.grid(True)
    plt.legend()
    plt.show()


print("--- Script Finished ---")