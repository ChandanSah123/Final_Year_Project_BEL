import os
import sys
import textwrap
import matplotlib.pyplot as plt
import config  # Imports your config.py

# ==============================================================================
# PSS/E INITIALIZATION
# ==============================================================================
try:
    import psse3603 # type: ignore
    import psspy    # type: ignore
    import dyntools # type: ignore
except ImportError:
    # Manual Fallback
    sys.path.append(r"C:\Program Files\PTI\PSSE35\PSSPY37")
    sys.path.append(r"C:\Program Files\PTI\PSSE35\PSSBIN")
    import psspy
    import dyntools

_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(1000)

# ==============================================================================
# SETUP
# ==============================================================================
raw_file = config.RAW_FILE
dyr_file = config.DYR_FILE
out_file = config.OUT_FILE.replace('.out', '_BusFault.out')

print("--- Simulation: Bus Fault (No Trip) ---")

# 1. Load Case
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")
psspy.fnsl([0,0,0,1,1,0,99,0])

# 2. Define Channels (Using "All" flag = 1)
psspy.delete_all_plot_channels()
# Record All Angles
psspy.chsb(0, 1, [-1, -1, -1, 1, 1, 0])
# Record All Speeds
psspy.chsb(0, 1, [-1, -1, -1, 1, 7, 0])
# Record All Active Powers
psspy.chsb(0, 1, [-1, -1, -1, 1, 2, 0])

# 3. Convert Network
psspy.cong(0)
psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])

# 4. Simulation Loop
psspy.strt(0, out_file)
psspy.run(0, config.T_FAULT_START, 0, 1, 0)

print(f"Applying Fault at Bus {config.FAULT_BUS}...")
psspy.dist_bus_fault(config.FAULT_BUS, 1, 0.0, [0.0, -0.2E+10])

t_clear = config.T_FAULT_START + config.T_CLEAR
psspy.run(0, t_clear, 0, 1, 0)

print("Clearing Fault (No Line Trip)...")
psspy.dist_clear_fault(1)

psspy.run(0, config.T_END, 0, 1, 0)
print("Simulation Complete.")

# ==============================================================================
# PLOTTING
# ==============================================================================
chnf_obj = dyntools.CHNF(out_file)
short_title, chanid, chandata = chnf_obj.get_data()
time = chandata['time']

fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
colors = ['r', 'g', 'b']

# Map Channels: PSS/E records them sequentially (Ang1,Ang2,Ang3, Spd1,Spd2,Spd3...)
# Assuming 3 Generators:
num_gens = len(config.GEN_BUSES)
for i, bus in enumerate(config.GEN_BUSES):
    c = colors[i % 3]
    lbl = f"Gen {bus}"
    
    # Calculate IDs
    id_ang = 1 + i
    id_spd = 1 + num_gens + i
    id_pe  = 1 + (2 * num_gens) + i
    
    if id_ang in chandata: axes[0].plot(time, chandata[id_ang], label=lbl, color=c)
    if id_spd in chandata: axes[1].plot(time, chandata[id_spd], label=lbl, color=c)
    if id_pe in chandata:  axes[2].plot(time, chandata[id_pe], label=lbl, color=c)

axes[0].set_title("Rotor Angle (Bus Fault)")
axes[0].legend()
axes[0].grid(True)
axes[1].set_title("Speed Deviation")
axes[1].grid(True)
axes[2].set_title("Active Power")
axes[2].grid(True)
plt.tight_layout()
plt.show()