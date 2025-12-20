import psse3603
import psspy
import dyntools
import matplotlib.pyplot as plt
import pandas as pd
import scipy.io
import numpy as np
import os
import re

# ==============================================================================
# 1. SETUP & CONFIGURATION
# ==============================================================================
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\WSCC_Programs_MY"
result_dir=r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\WSCC_Programs_MY\Result_Directory"
raw_file = os.path.join(work_dir, "IEEE9bus.raw")
dyr_file = os.path.join(work_dir, "ieee9bus.dyr")
out_file = os.path.join(result_dir, "IEEE9.out")
excel_file = os.path.join(result_dir, "IEEE9_Results.xlsx")
mat_file = os.path.join(result_dir, "IEEE9_Results.mat")

# Initialize PSS/E
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
print("--- Starting Simulation for IEEE 9-Bus ---")

# Simulation Parameters
fault_bus = 7
from_bus=7
to_bus=5
line_id='1'
fault_time = 1
fault_duration=0.15
clear_time = fault_time+fault_duration
end_time = 4

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
    plt.title("Rotor Angles")
    plt.ylabel("Angle (Deg)")
    plt.grid(True)
    plt.legend()
    plt.show(block=False)

    # Plot 2: Mechanical Power (Channels 4, 5, 6)
    plt.figure(figsize=(10, 6))
    for ch in range(4, 7):
        if ch in chandata:
            clean_label = str(chanid[ch]).strip()
            plt.plot(time_data, chandata[ch], label=f"Ch{ch}: {clean_label}")
    plt.title("Mechanical Power (Pmech)")
    plt.ylabel("Power (pu)")
    plt.grid(True)
    plt.legend()
    plt.show(block=False)

    # Plot 3: Electrical Power (Channels 7, 8, 9)
    plt.figure(figsize=(10, 6))
    for ch in range(7, 10):
        if ch in chandata:
            clean_label = str(chanid[ch]).strip()
            plt.plot(time_data, chandata[ch], label=f"Ch{ch}: {clean_label}")
    plt.title("Electrical Power (Pelec)")
    plt.ylabel("Power (pu)")
    plt.grid(True)
    plt.legend()
    plt.show(block=False)

    # Plot 4: Generator Speed (Channels 10, 11, 12)
    plt.figure(figsize=(10, 6))
    for ch in range(10, 13):
        if ch in chandata:
            clean_label = str(chanid[ch]).strip()
            plt.plot(time_data, chandata[ch], label=f"Ch{ch}: {clean_label}")
    plt.title("Generator Speed")
    plt.ylabel("Speed (pu)")
    plt.grid(True)
    plt.legend()
    plt.show() # Block here to keep windows open

print("--- Script Finished ---")