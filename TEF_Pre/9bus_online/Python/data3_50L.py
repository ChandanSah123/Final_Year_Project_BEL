import os
import sys
import textwrap
import scipy.io
import numpy as np
import config

try:
    import psse3603; import psspy; import dyntools
except ImportError:
    sys.path.append(r"C:\Program Files\PTI\PSSE35\PSSPY37")
    sys.path.append(r"C:\Program Files\PTI\PSSE35\PSSBIN")
    import psspy; import dyntools

_i = psspy.getdefaultint(); _f = psspy.getdefaultreal()
psspy.psseinit(1000)

# ==============================================================================
# 1. CONFIGURATION: 50% LINE FAULT (TRIP)
# ==============================================================================
raw_file = config.RAW_FILE
dyr_file = config.DYR_FILE
mat_file = config.MAT_FILE_DATA
out_file = config.OUT_FILE.replace('.out', '_Line50.out')

from_bus = config.TRIP_LINE_FROM
to_bus   = config.TRIP_LINE_TO
ckt_id   = config.LINE_ID
t_fault  = config.T_FAULT_START
t_clear  = config.T_FAULT_START + config.T_CLEAR
t_end    = config.T_END

print("--- Starting Simulation: 50% Line Fault (Trip) ---")

# ==============================================================================
# 2. SIMULATION
# ==============================================================================
psspy.read(0, raw_file)
psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")
psspy.fnsl([0,0,0,1,1,0,99,0])

# Define Channels
psspy.delete_all_plot_channels()
psspy.chsb(0, 1, [-1, -1, -1, 1, 1, 0]) # Angle
psspy.chsb(0, 1, [-1, -1, -1, 1, 7, 0]) # Speed
psspy.chsb(0, 1, [-1, -1, -1, 1, 6, 0]) # Pm
psspy.chsb(0, 1, [-1, -1, -1, 1, 2, 0]) # Pe
psspy.chsb(0, 1, [-1, -1, -1, 1, 4, 0]) # Q

# Convert
psspy.cong(0)
psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])
psspy.fact(); psspy.tysl(0)

# Run
psspy.strt(0, out_file)
psspy.run(0, t_fault, 0, 1, 0)

print(f"Applying 50% Fault Proxy at Bus {from_bus}...")
# NOTE: Using standard bus fault as proxy for 50% line fault behavior
psspy.dist_bus_fault(from_bus, 1, 0.0, [0.0, -0.2E+10])

psspy.run(0, t_clear, 0, 1, 0)

print(f"Clearing Fault by Tripping Line {from_bus}-{to_bus}...")
psspy.dist_branch_trip(from_bus, to_bus, ckt_id)
psspy.dist_clear_fault(1)

psspy.run(0, t_end, 0, 1, 0)

# ==============================================================================
# 3. EXTRACTION
# ==============================================================================
print("Extracting Data...")
chnf_obj = dyntools.CHNF(out_file)
short_title, chanid, chandata = chnf_obj.get_data()
t_data = chandata['time']
sorted_keys = sorted([k for k in chanid.keys() if isinstance(k, int)])
data_columns = []

for key in sorted_keys:
    y_data = chandata[key]
    data_columns.append(y_data[:len(t_data)])

if data_columns:
    time_col = t_data[:len(data_columns[0])]
    data_columns.append(time_col)
    final_matrix = np.column_stack(data_columns)
    scipy.io.savemat(mat_file, {"data": final_matrix})
    print(f"SUCCESS: {mat_file} saved.")
    print(f"Matrix Shape: {final_matrix.shape}")