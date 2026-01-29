# config.py
import os

# ==============================================================================
# 1. PATHS & DIRECTORIES
# ==============================================================================
BASE_DIR = os.path.dirname(__file__)
WORK_DIR = BASE_DIR
RESULT_DIR = os.path.join(os.path.dirname(BASE_DIR), "Matlab script")


# Auto-generate full file paths
RAW_FILE = os.path.join(WORK_DIR, "ieee39bus1.raw")
DYR_FILE = os.path.join(WORK_DIR, "ieee39buscls.dyr")
OUT_FILE = os.path.join(RESULT_DIR, "IEEE39_Combined.out")
MAT_FILE_DATA = os.path.join(RESULT_DIR, "data1.mat")
MAT_FILE_YBUS = os.path.join(RESULT_DIR, "Y_all.mat")
MAT_FILE_CCT  = os.path.join(RESULT_DIR, "CCT_TimeDomain.mat")
TXT_FILE_YBUS = os.path.join(RESULT_DIR, "Ybus_Export.txt")
TIME_FILE = os.path.join(RESULT_DIR, "Tcl.mat")

# ==============================================================================
# 2. FAULT PARAMETERS (Change these to update ALL scripts)
# ==============================================================================
FAULT_BUS = 1
TRIP_LINE_FROM = 1
TRIP_LINE_TO = 2
LINE_ID = '1 '

# Timing
T_FAULT_START = 1.0
# Single centralized clearing time used by all scripts
T_CLEAR = 1
T_END = T_FAULT_START + T_CLEAR+3
# ==============================================================================
# 3. SYSTEM DATA (For Y-Bus Calculation)
# ==============================================================================
GEN_BUSES = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39] 
XD_PRIME = [0.025, 0.05, 0.045, 0.035, 0.089, 0.04, 0.044, 0.045, 0.045, 0.004]

# Hardcoded Loads (Exact match to MATLAB/System)
LOAD_ADM = {
  
                1:	0,
                2:	0,
                3:	3.034-0.0226j,
                4:	4.9612-1.8257j,
                5:	0,
                6:	0,
                7:	2.3521-0.8451j,
                8:	5.262-1.7742j,
                9:	0,
                10:	0,
                11:	0,
                12:	0,
                13:	0,
                14:	0,
                15:	3.1037-1.4839j,
                16:	3.0903-0.3034j,
                17:	0,
                18:	1.4867-0.2823j,
                19:	0,
                20:	6.392-1.0484j,
                21:	2.5737-1.0802j,
                22:	0,
                23:	2.2673-0.775j,
                24:	2.8681+0.855j,
                25:	2.0027-0.422j,
                26:	1.2557-0.1536j,
                27:	2.6095-0.7011j,
                28:	1.8681-0.2503j,
                29:	2.5719-0.244j,
                31:	0.0838-0.0419j,
                39:	10.4063-2.3565j

}