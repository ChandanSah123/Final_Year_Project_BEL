# config.py
import os

# ==============================================================================
# 1. PATHS & DIRECTORIES
# ==============================================================================
BASE_DIR = os.path.dirname(__file__)
WORK_DIR = BASE_DIR
RESULT_DIR = os.path.join(os.path.dirname(BASE_DIR), "Matlab")


# Auto-generate full file paths
RAW_FILE = os.path.join(WORK_DIR, "IEEE14bus_raw.raw")
DYR_FILE = os.path.join(WORK_DIR, "ieee14.dyr")
OUT_FILE = os.path.join(RESULT_DIR, "IEEE14.out")
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
TRIP_LINE_TO = 5
LINE_ID = '1 '

# Timing
T_FAULT_START = 1
# Single centralized clearing time used by all scripts
T_CLEAR = 0.12
T_END = T_FAULT_START + T_CLEAR+3
# ==============================================================================
# 3. SYSTEM DATA (For Y-Bus Calculation)
# ==============================================================================
GEN_BUSES = [1, 2, 3, 6, 8]
#XD_PRIME = [0.30, 0.19, 0.185, 0.232, 0.232]
XD_PRIME = [0.230, 0.13, 0.13, 0.12, 0.12]

# Hardcoded Loads (Exact match to MATLAB/System)
LOAD_ADM = {
    
2:0.1987-0.1163j,
3:0.9234-0.1863j,
4:0.4607,
5:0.072-0.0152j,
6:0.1029-0.0689j,
9:0.2863-0.1611j,
10:0.0878-0.0566j,
11:0.0334-0.0172j,
12:0.0576-0.0151j,
13:0.1282-0.0551j,
14:0.1432-0.0481j


}