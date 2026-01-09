import os

# Base directories (auto-detect relative to this file)
BASE_DIR = os.path.dirname(__file__)
WORK_DIR = BASE_DIR
RESULT_DIR = os.path.join(os.path.dirname(BASE_DIR), "Matlab script")

# Common file names
DYR_FILENAME = "ieee39buscls.dyr"
DYR_FILE = os.path.join(WORK_DIR, DYR_FILENAME)

# PSS/E initialization settings
# Number of buses for psseinit. Use a value large enough for the system.
PSSE_INIT = 50

# Dynamics solver timestep
DYNAMICS_STEP = 0.001

# Fault / scenario defaults (centralized)
FAULT_BUS = 1
TRIP_LINE_FROM = 1
TRIP_LINE_TO = 2
CKT_ID = '1 '
FAULT_TIME = 1.0
# Default clearing time used by some scripts (can be overridden per-run)
DEFAULT_CLEARING = 0.32
# Generic helpers
def out_path(filename):
    return os.path.join(RESULT_DIR, filename)

def raw_path(filename):
    return os.path.join(WORK_DIR, filename)
