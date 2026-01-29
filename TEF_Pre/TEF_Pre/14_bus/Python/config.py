# config.py
import os

# ==============================================================================
# 1. PATHS & DIRECTORIES
# ==============================================================================
BASE_DIR = os.path.dirname(__file__)
WORK_DIR = BASE_DIR
RESULT_DIR = os.path.join(os.path.dirname(BASE_DIR), "Matlab")


# Auto-generate full file paths
RAW_FILE = os.path.join(WORK_DIR, "IEEE9bus.raw")
DYR_FILE = os.path.join(WORK_DIR, "ieee9bus.dyr")
OUT_FILE = os.path.join(RESULT_DIR, "IEEE9_Combined.out")
MAT_FILE_DATA = os.path.join(RESULT_DIR, "data1.mat")
MAT_FILE_YBUS = os.path.join(RESULT_DIR, "Y_all.mat")
MAT_FILE_CCT  = os.path.join(RESULT_DIR, "CCT_TimeDomain.mat")
TXT_FILE_YBUS = os.path.join(RESULT_DIR, "Ybus_Export.txt")
TIME_FILE = os.path.join(RESULT_DIR, "Tcl.mat")

# ==============================================================================
# 2. FAULT PARAMETERS (Change these to update ALL scripts)
# ==============================================================================
FAULT_BUS = 7
TRIP_LINE_FROM = 7
TRIP_LINE_TO = 5
LINE_ID = '1 '

# Timing
T_FAULT_START = 1.0
# Single centralized clearing time used by all scripts
T_CLEAR = 0.163
T_END = T_FAULT_START + T_CLEAR+3
# ==============================================================================
# 3. SYSTEM DATA (For Y-Bus Calculation)
# ==============================================================================
GEN_BUSES = [1, 2, 3]
XD_PRIME = [0.0608, 0.1198, 0.1813]

# Hardcoded Loads (Exact match to MATLAB/System)
LOAD_ADM = {
    5: 1.26 - 0.504j,
    6: 0.877 - 0.292j,
    8: 0.969 - 0.339j
}