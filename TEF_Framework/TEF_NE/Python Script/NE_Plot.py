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
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\NE"
raw_file = os.path.join(work_dir, "ieee39bus1.raw")
dyr_file = os.path.join(work_dir, "ieee39buscls1.dyr")
out_file = os.path.join(work_dir, "NE.out")

# Initialize PSS/E
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
print("--- Starting Simulation ---")

# Simulation Parameters
tf = [0.87, 0.328, 0.306, 0.311, 0.298, 0.29, 0.35, 0.342, 1, 0.293, 0.309, 1, 0.318, 0.314, 0.276, 0.203, 0.242, 0.293, 0.245, 0.284, 0.238, 0.21, 0.197, 0.256, 0.27, 0.213, 0.292, 0.215, 0.183, 1, 0.323, 0.294, 0.272, 0.275, 0.262, 0.189, 0.292, 0.176]
b = 22
fault_time = 1
clear_time = 1.3
end_time = 10

# Load Case
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

# Solver Settings
psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
psspy.fnsl([0,0,0,1,1,0,99,0])

# Define Channels (What to record)
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,2,0])
psspy.chsb(0,1,[-1,-1,-1,1,3,0])
psspy.chsb(0,1,[-1,-1,-1,1,7,0])

# Convert Network
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])

# Run Simulation
psspy.strt(0, out_file)
psspy.run(0, fault_time, 0, 1, 0)

print(f"Applying Fault at Bus {b}...")
psspy.dist_bus_fault(b, 1, 0.0, [0.0, -0.2E+10])

psspy.run(0, clear_time, 0, 1, 0)

print("Clearing Fault...")
psspy.dist_clear_fault(1)

psspy.run(0, end_time, 0, 1, 0)
print("Simulation Complete.")

# ==============================================================================
# 3. DATA EXTRACTION
# ==============================================================================
chnf_obj = dyntools.CHNF(out_file)
short_title, chanid, chandata = chnf_obj.get_data()
time_data = chandata['time']

# ==============================================================================
# 4. PLOTTING
# ==============================================================================
print("Generating Plots...")

# Plot 1: Rotor Angles (1-10)
plt.figure(figsize=(10, 6))
for ch in range(1, 11):
    if ch in chandata:
        clean_label = str(chanid[ch]).strip()
        plt.plot(time_data, chandata[ch], label=clean_label)
plt.title("Rotor Angles")
plt.xlabel("Time (s)")
plt.ylabel("Angle (Deg)")
plt.grid(True)
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
plt.tight_layout()
plt.show(block=False)

# Plot 2: Speed (31-40)
plt.figure(figsize=(10, 6))
for ch in range(31, 41):
    if ch in chandata:
        clean_label = str(chanid[ch]).strip()
        plt.plot(time_data, chandata[ch], label=clean_label)
plt.title("Generator Speed")
plt.xlabel("Time (s)")
plt.ylabel("Speed (pu)")
plt.grid(True)
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
plt.tight_layout()
plt.show()

