import os
import sys
import textwrap  # CRITICAL FIX: Required by dyntools in some PSS/E versions
import matplotlib.pyplot as plt
import config  # Imports your config.py

# ==============================================================================
# 1. SCENARIO CONFIGURATION
# ==============================================================================
# 1 = Bus Fault (No Trip), 2 = Line Fault (Trip), 3 = 50% Line Fault (Trip)
FAULT_SCENARIO = 2  

# ==============================================================================
# 2. PSS/E INITIALIZATION
# ==============================================================================
try:
    import psse3603 # type: ignore
    import psspy    # type: ignore
    import dyntools # type: ignore
    pss_available = True
except ImportError:
    # Fallback paths
    sys.path.append(r"C:\Program Files\PTI\PSSE35\PSSPY37")
    sys.path.append(r"C:\Program Files\PTI\PSSE35\PSSBIN")
    try:
        import psspy
        import dyntools
        pss_available = True
    except ImportError:
        print("CRITICAL: PSS/E libraries not found.")
        sys.exit()

_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(1000)

# ==============================================================================
# 3. SETUP
# ==============================================================================
raw_file = config.RAW_FILE
dyr_file = config.DYR_FILE
out_file = config.OUT_FILE.replace('.out', f'_Scen{FAULT_SCENARIO}.out')

# Parameters
fault_bus = config.FAULT_BUS
from_bus  = config.TRIP_LINE_FROM
to_bus    = config.TRIP_LINE_TO
ckt_id    = config.LINE_ID
gen_id    = '1 ' 

t_fault = config.T_FAULT_START
t_clear = config.T_FAULT_START + config.T_CLEAR
t_end   = config.T_END

scen_name = {1: "Bus Fault (No Trip)", 2: "Line Fault (Trip)", 3: "50% Fault (Trip)"}
print(f"--- Simulation: {scen_name.get(FAULT_SCENARIO)} ---")

# ==============================================================================
# 4. SIMULATION
# ==============================================================================
# A. Load & Power Flow
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")
psspy.fnsl([0,0,0,1,1,0,99,0])

# B. Convert Network (Norton)
psspy.cong(0)
psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])

# C. DEFINE CHANNELS (CORRECTED)
psspy.delete_all_plot_channels()

print("Defining Channels...")
for bus in config.GEN_BUSES:
    # FIX: Status array must be [1, -1, -1, 1, 1, 0] or similar standard format.
    # The API 'machine_array_channel' is tricky. 
    # Let's use the safer 'chsb' (Channel Select Bus) which you used originally but logic-ed wrong.
    # We will loop explicitly.
    
    # 1 = ANGLE
    # chsb args: (sid, all, status[6])
    # status: [-1, -1, -1, 1, 1, 0] selects Angle (1) for machine at 'bus'
    ierr = psspy.chsb(0, 0, [-1, -1, -1, 1, 1, 0], bus, gen_id)
    if ierr != 0: print(f"Error defining Angle for Bus {bus}")

    # 7 = SPEED
    # status: [-1, -1, -1, 1, 7, 0] selects Speed (7)
    ierr = psspy.chsb(0, 0, [-1, -1, -1, 1, 7, 0], bus, gen_id)
    if ierr != 0: print(f"Error defining Speed for Bus {bus}")

    # 2 = PELEC
    # status: [-1, -1, -1, 1, 2, 0] selects Pe (2)
    ierr = psspy.chsb(0, 0, [-1, -1, -1, 1, 2, 0], bus, gen_id)
    if ierr != 0: print(f"Error defining Pe for Bus {bus}")

# D. Initialize State
ierr_strt = psspy.strt(0, out_file)
if ierr_strt > 0:
    print(f"CRITICAL: Simulation failed to start (Error {ierr_strt})")
    sys.exit()

# E. Time Stepping
psspy.run(0, t_fault, 0, 1, 0)

print(f"Applying Fault at t={t_fault}")
if FAULT_SCENARIO == 1:
    psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])
elif FAULT_SCENARIO == 2:
    psspy.dist_bus_fault(from_bus, 1, 0.0, [0.0, -0.2E+10]) 
elif FAULT_SCENARIO == 3:
    psspy.dist_bus_fault(from_bus, 1, 0.0, [0.0, -0.2E+10])

psspy.run(0, t_clear, 0, 1, 0)

print(f"Clearing Fault at t={t_clear}")
if FAULT_SCENARIO == 1:
    psspy.dist_clear_fault(1)
else:
    psspy.dist_branch_trip(from_bus, to_bus, ckt_id)
    psspy.dist_clear_fault(1)

psspy.run(0, t_end, 0, 1, 0)
print("Simulation Finished. Processing Data...")

# ==============================================================================
# 5. DATA EXTRACTION & PLOTTING
# ==============================================================================
if not os.path.exists(out_file):
    print("Error: .out file was not created.")
    sys.exit()

# Setup dyntools
try:
    chnf_obj = dyntools.CHNF(out_file)
    short_title, chanid, chandata = chnf_obj.get_data()
except Exception as e:
    print(f"Error reading .out file: {e}")
    sys.exit()

# Debug: Check if data exists
if 'time' not in chandata or len(chandata['time']) == 0:
    print("ERROR: No data recorded in .out file (Simulation might have crashed immediately).")
    sys.exit()
else:
    print(f"Data Loaded: {len(chandata['time'])} time points.")
    # Show first few keys to help debug channel mapping
    keys = list(chandata.keys())
    print(f"Channel Keys present: {keys[:15]}...") 

time = chandata['time']
fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

# Colors for generators
colors = ['r', 'g', 'b']

# Note on Channel IDs with 'chsb':
# IDs are assigned sequentially as we called 'chsb'.
# We called Angle, Speed, Pe for Bus 1, then Bus 2, then Bus 3.
# So IDs should be:
# 1: Ang(1), 2: Spd(1), 3: Pe(1)
# 4: Ang(2), 5: Spd(2), 6: Pe(2)
# 7: Ang(3), 8: Spd(3), 9: Pe(3)

for i, bus in enumerate(config.GEN_BUSES):
    # Calculate Channel ID (1-based index)
    base_id = (i * 3) + 1 
    id_ang = base_id
    id_spd = base_id + 1
    id_pe  = base_id + 2
    
    label = f"Gen {bus}"
    c = colors[i % len(colors)]
    
    # Plot Angle
    if id_ang in chandata:
        axes[0].plot(time, chandata[id_ang], label=label, color=c)
    
    # Plot Speed
    if id_spd in chandata:
        axes[1].plot(time, chandata[id_spd], label=label, color=c)
        
    # Plot Power
    if id_pe in chandata:
        axes[2].plot(time, chandata[id_pe], label=label, color=c)

# Formatting
axes[0].set_title(f"Rotor Angle\n{scen_name.get(FAULT_SCENARIO)}")
axes[0].set_ylabel("Angle (deg)")
axes[0].grid(True)
axes[0].legend(loc="upper left", fontsize='small')

axes[1].set_title("Speed Deviation")
axes[1].set_ylabel("Speed (pu)")
axes[1].grid(True)

axes[2].set_title("Active Power")
axes[2].set_ylabel("Power (MW/pu)")
axes[2].set_xlabel("Time (s)")
axes[2].grid(True)

plt.tight_layout()
print("Displaying Plot...")
plt.show()