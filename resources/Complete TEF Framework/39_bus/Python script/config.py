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
FAULT_BUS = 29
TRIP_LINE_FROM = 29
TRIP_LINE_TO =26
LINE_ID = '1 '

# Timing
T_FAULT_START = 1.0
# Single centralized clearing time used by all scripts
T_CLEAR = 0.184
T_END = 5
# ==============================================================================
# 3. SYSTEM DATA (For Y-Bus Calculation)
# ==============================================================================
GEN_BUSES = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39] 
XD_PRIME = [0.025, 0.05, 0.045, 0.035, 0.089, 0.04, 0.044, 0.045, 0.045, 0.004]

# Hardcoded Loads (Exact match to MATLAB/System)
LOAD_ADM = {
  
                3:3.0593-0.0228j,
                4:4.9995-1.8398j,
                7:2.3636-0.8492j,
                8:5.2871-1.7826j,
                12:0.0754-0.885j,
                15:3.1574-1.5096j,
                16:3.1529-0.3067j,
                18:1.5054-0.2858j,
                20:6.419-1.0528j,
                21:2.6431-1.1093j,
                23:2.3978-0.8196j,
                24:2.9336+0.8746j,
                25:2.0098-0.4235j,
                26:1.2648-0.1547j,
                27:2.638-0.7088j,
                28:1.875-0.2512j,
                29:2.5782-0.2446j,
                31:0.0954-0.0477j,
                39:10.4063-2.3565j

}