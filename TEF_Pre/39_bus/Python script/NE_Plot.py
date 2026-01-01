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
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_Pre\39_bus\Python Script"
result_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_Pre\39_bus\Matlab script"
raw_file = os.path.join(work_dir, "IEEE39bus1.raw")
dyr_file = os.path.join(work_dir, "ieee39buscls.dyr")
out_file = os.path.join(result_dir, "IEEE39.out")

# Initialize PSS/E
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
print("--- Starting Simulation for IEEE 39-Bus ---")

# Simulation Parameters

fault_bus = 36
from_bus = 36
to_bus = 23
line_id = '1 '
fault_time = 1
t_cl=0.0
clear_time = fault_time+t_cl
end_time = 5

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
# 3. PLOTTING (Distinct Variations for All Generators)
# ==============================================================================

if not os.path.exists(out_file):
    print(f"CRITICAL ERROR: Output file not found at {out_file}")
else:
    # 1. Load Data
    chnf_obj = dyntools.CHNF(out_file)
    short_title, chanid, chandata = chnf_obj.get_data()
    time_data = list(chandata['time'])

    print("Generating Distinct Generator Plots...")

    # 2. Setup Figure
    plt.figure(figsize=(12, 8)) # Use a large figure size for clarity

    # 3. Define Styles to Cycle Through
    # We mix 4 line styles with standard colors to ensure 10 distinct combos
    linestyles = ['-', '--', '-.', ':'] 
    
    # Get the default matplotlib color cycle (usually 10 distinct colors)
    # This ensures Gen 1 is Blue, Gen 2 is Orange, etc.
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color'] 

    # 4. Loop Through All Generator Channels (1-10)
    for i, ch in enumerate(range(1, 11)):
        if ch in chandata:
            # Get Data
            y_data = chandata[ch]
            clean_label = str(chanid[ch]).strip()
            
            # --- ASSIGN UNIQUE VARIATION ---
            # Cycle through colors (0-9)
            color = colors[i % len(colors)]
            
            # Cycle through line styles (Solid, Dashed, Dash-Dot, Dotted)
            # shifting style every 3 generators prevents "Blue Solid" and "Blue Dashed" looking too similar
            style = linestyles[i % len(linestyles)]
            
            # PLOT (Standard width, no fading)
            plt.plot(time_data, y_data, 
                     label=f"{clean_label}", 
                     color=color,
                     linestyle=style,
                     linewidth=2.0) # Slightly thicker for visibility

    # 5. Formatting
    plt.title(f"Generator Rotor Angles (t_cl = {t_cl}s)", fontsize=14)
    plt.ylabel("Angle (Deg)", fontsize=12)
    plt.xlabel("Time (s)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Move Legend outside the plot so it doesn't cover the data
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title="Generators")
    
    plt.tight_layout() # Adjusts margins to fit the legend
    plt.show()

print("--- Script Finished ---")