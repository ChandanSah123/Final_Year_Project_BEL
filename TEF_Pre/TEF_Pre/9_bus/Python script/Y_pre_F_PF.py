import os
import re
import numpy as np
import scipy.io as sio
import sys

# Attempt to import PSS/E modules (safe)
try:
    import psse3603  # type: ignore
    import psspy     # type: ignore
    import dyntools  # type: ignore
    pss_available = True
except Exception:
    pss_available = False

import config

# ==============================================================================
# 1. SCENARIO CONFIGURATION (IEEE 9-BUS)
# ==============================================================================
work_dir = config.WORK_DIR
result_dir = config.RESULT_DIR

raw_file = config.RAW_FILE
txt_file = config.TXT_FILE_YBUS
mat_file = config.MAT_FILE_YBUS

# --- SCENARIO SETTINGS ---
# To replicate your specific case: "Fault at Bus 7, Trip Line 7-5"
FAULT_TYPE = 'BUS'          # 'BUS' or 'LINE'
CLEARING_ACTION = 'TRIP'    # 'TRIP' or 'NONE'

# IF FAULT_TYPE = 'BUS':
FAULT_BUS_ID = config.FAULT_BUS

# IF FAULT_TYPE = 'LINE' or CLEARING_ACTION = 'TRIP':
FAULT_LINE_FROM = config.TRIP_LINE_FROM
FAULT_LINE_TO = config.TRIP_LINE_TO
FAULT_LINE_ID = config.LINE_ID
FAULT_RATIO = 0.5       # (Ignored for Bus Faults)

# --- SYSTEM PARAMETERS (IEEE 9-Bus) ---
gen_buses = [1, 2, 3]
xd_prime = [0.0608, 0.1198, 0.1813]

# Hardcoded Loads (Constant Impedance)
load_adm = {
    5: 1.26 - 0.504j,
    6: 0.877 - 0.292j,
    8: 0.969 - 0.339j
}

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

def init_psse_and_read_case():
    try:
        _ = psspy.psseinit(200)
        ierr = psspy.read(0, raw_file)
        if ierr != 0:
            raise RuntimeError(f"psspy.read returned error {ierr}")
    except Exception:
        raise

def get_line_params_from_psse(bus_from, bus_to, ckt_id='1 '):
    """ Extracts Y_series and B_charging/magnetizing from PSS/E. """
    z_complex = 0j
    b_total_real = 0.0
    found = False

    # Try Line
    ierr, z_val = psspy.brndt2(bus_from, bus_to, ckt_id, "RX")
    if ierr == 0:
        ierr_c, b_val = psspy.brndat(bus_from, bus_to, ckt_id, "CHARG")
        if ierr_c == 0:
            z_complex = z_val
            b_total_real = b_val
            found = True
            print(f"Detected Line {bus_from}-{bus_to}")

    # Try Transformer
    if not found:
        ierr, z_val = psspy.xfrdat(bus_from, bus_to, ckt_id, "RX")
        if ierr == 0:
            ierr_m, mag_imag = psspy.xfrdat(bus_from, bus_to, ckt_id, "MAG2")
            z_complex = z_val
            b_total_real = mag_imag if ierr_m == 0 else 0.0
            found = True
            print(f"Detected Transformer {bus_from}-{bus_to}")

    if not found:
        # Check reverse just in case
        raise ValueError(f"Could not find branch {bus_from}-{bus_to} ID {ckt_id}")

    if z_complex == 0j: raise ValueError(f"Zero Impedance for {bus_from}-{bus_to}")
    
    Y_line = 1.0 / z_complex
    B_line_half = 1j * (b_total_real / 2.0)
    return Y_line, B_line_half

def solve_and_export_ybus():
    """ Solves power flow and exports the physical Ybus (Lines/Trafo/Shunts). """
    psspy.fnsl([0,0,0,1,1,0,99,0])
    psspy.cong(0) # Norton equivalent for generators
    # No conl() -> We add loads manually
    try:
        psspy.ordr(0)
        psspy.fact()
        psspy.tysl(0)
    except: pass
    
    ierr = psspy.output_y_matrix(0, 1, 0, 0, txt_file)
    if isinstance(ierr, tuple): ierr = ierr[0]
    if ierr != 0: raise RuntimeError(f"Ybus export failed: {ierr}")

def parse_ybus_text_to_numpy():
    """ Reads PSS/E Ybus text file into Numpy Matrix. """
    Yreal, Yimag = {}, {}
    max_bus = 0
    float_re = re.compile(r"[-+]?\d*\.\d+|\d+")
    
    with open(txt_file, "r") as f:
        for line in f:
            parts = float_re.findall(line)
            if len(parts) >= 4:
                try:
                    i, j = int(parts[0]), int(parts[1])
                    real, imag = float(parts[2]), float(parts[3])
                    max_bus = max(max_bus, i, j)
                    Yreal[(i, j)] = real
                    Yimag[(i, j)] = imag
                except: continue
    
    Y = np.zeros((max_bus, max_bus), dtype=np.complex128)
    for (i, j), rv in Yreal.items():
        iv = Yimag.get((i, j), 0.0)
        Y[i-1, j-1] = rv + 1j*iv
        Y[j-1, i-1] = rv + 1j*iv
    return Y

def add_loads_to_Ybus(Ybus, load_admittance):
    """ Adds load admittances to the physical matrix diagonals """
    Y_loaded = Ybus.copy()
    for bus_num, yload in load_admittance.items():
        if bus_num <= Y_loaded.shape[0]:
            idx = bus_num - 1
            Y_loaded[idx, idx] += yload
    return Y_loaded

# --- MATRIX MODIFICATION LOGIC ---

def build_faulted_bus_matrix(Y_full, fault_bus_num):
    """ SCENARIO 1: Fault at a Bus (Ground the bus -> Remove Row/Col) """
    idx = fault_bus_num - 1
    Yf = np.delete(Y_full, idx, axis=0)
    Yf = np.delete(Yf, idx, axis=1)
    return Yf

def build_faulted_line_matrix(Y_full, from_bus, to_bus, y_series, ratio=0.5):
    """ SCENARIO 2: Fault on a Line (Split line and ground middle) """
    idx_a = from_bus - 1
    idx_b = to_bus - 1
    Yf = Y_full.copy()
    
    # 1. Remove the original line A-B from the matrix
    # Subtract y_series from diagonals A and B
    Yf[idx_a, idx_a] -= y_series
    Yf[idx_b, idx_b] -= y_series
    
    # Remove from off-diagonal (Add y_series because Yij = -y_series)
    Yf[idx_a, idx_b] += y_series 
    Yf[idx_b, idx_a] += y_series

    # 2. Add Shunt Admittances for the Fault (Pi-Split)
    y_shunt_a = 1.0 / ( (1.0/y_series) * ratio )
    y_shunt_b = 1.0 / ( (1.0/y_series) * (1.0 - ratio) )
    
    Yf[idx_a, idx_a] += y_shunt_a
    Yf[idx_b, idx_b] += y_shunt_b
    
    return Yf

def build_tripped_line_matrix(Y_full, from_bus, to_bus, y_line, b_charging):
    """ Post-Fault: Line Tripped (Remove Line) """
    idx_a = from_bus - 1
    idx_b = to_bus - 1
    Yp = Y_full.copy()
    
    # Total admittance to subtract (Series + Charging)
    y_total = y_line + b_charging
    
    # Subtract from diagonals
    Yp[idx_a, idx_a] -= y_total
    Yp[idx_b, idx_b] -= y_total
    
    # Remove from off-diagonals (Add y_line because Yij = -y_line)
    Yp[idx_a, idx_b] += y_line
    Yp[idx_b, idx_a] += y_line
    
    return Yp

def kron_reduce_internal(Y_physical, gen_buses_list, xd_prime_list, load_admittance, is_bus_fault_active=False, fault_bus_idx=None):
    """
    Standard Kron Reduction.
    If is_bus_fault_active is True, it handles the shifted indices for buses > fault_bus.
    """
    K = len(gen_buses_list)
    N = Y_physical.shape[0]

    # A. Generator Admittances (ybar)
    ybar = np.zeros((K, K), dtype=complex)
    for i, xd in enumerate(xd_prime_list):
        ybar[i, i] = 1.0 / (1j * xd)

    YD = Y_physical.copy()

    # B. Add Gen & Load Admittances to Diagonals
    # Helper to find physical index
    def get_phys_idx(bus_id):
        if is_bus_fault_active:
            if bus_id == fault_bus_idx: return None # This bus doesn't exist anymore
            return bus_id - 1 if bus_id < fault_bus_idx else bus_id - 2
        else:
            return bus_id - 1

    gen_map = []
    for i, bus in enumerate(gen_buses_list):
        pidx = get_phys_idx(bus)
        if pidx is not None:
            gen_map.append((i, pidx))
            YD[pidx, pidx] += ybar[i, i]
            
    if load_admittance:
        for bus, yload in load_admittance.items():
            pidx = get_phys_idx(bus)
            if pidx is not None:
                YD[pidx, pidx] += yload

    # C. Partition and Solve
    YA = ybar.copy()
    YB = np.zeros((K, N), dtype=complex)
    for (k, pidx) in gen_map:
        YB[k, pidx] = -ybar[k, k]
    YC = YB.T

    try:
        # Yint = YA - YB * inv(YD) * YC
        # Use solve for stability: YD_inv_YC = YD \ YC
        Yint = YA - np.dot(YB, np.linalg.solve(YD, YC))
    except np.linalg.LinAlgError:
        print("Warning: Matrix singular during Kron Reduction.")
        Yint = np.zeros((K, K))
        
    return Yint

# ==============================================================================
# 3. MAIN EXECUTION
# ==============================================================================
if __name__ == "__main__":
    if not os.path.exists(raw_file):
        raise FileNotFoundError(f"Raw file not found: {raw_file}")
    
    print(f"--- SIMULATION SETUP (IEEE 9-BUS) ---")
    print(f"Fault Type:     {FAULT_TYPE}")
    if FAULT_TYPE == 'LINE':
        print(f"Location:       Line {FAULT_LINE_FROM}-{FAULT_LINE_TO} @ {FAULT_RATIO*100}%")
    else:
        print(f"Location:       Bus {FAULT_BUS_ID}")
    print(f"Clearing:       {CLEARING_ACTION}")
    print("-" * 30)

    # 1. Initialize PSS/E and Export Base Ybus
    init_psse_and_read_case()
    
    # Get Line Params if needed
    y_line_dyn = 0j
    b_line_dyn = 0j
    
    if FAULT_TYPE == 'LINE' or CLEARING_ACTION == 'TRIP':
        f_from = FAULT_LINE_FROM
        f_to   = FAULT_LINE_TO
        f_id   = FAULT_LINE_ID
        
        print(f"Extracting line parameters for {f_from}-{f_to}...")
        y_line_dyn, b_line_dyn = get_line_params_from_psse(f_from, f_to, f_id)

    solve_and_export_ybus()
    Y_pre_phys = parse_ybus_text_to_numpy() # 9x9 Matrix

    # ==========================================================================
    # STEP 1: PRE-FAULT MATRIX
    # ==========================================================================
    print("Calculating Yint_pre...")
    Yint_pre = kron_reduce_internal(Y_pre_phys, gen_buses, xd_prime, load_adm, is_bus_fault_active=False)

    # ==========================================================================
    # STEP 2: FAULT MATRIX (Y_fault)
    # ==========================================================================
    print("Calculating Yint_fault...")
    
    if FAULT_TYPE == 'BUS':
        # Remove Bus Row/Col (Matrix becomes 8x8)
        Y_phys_fault = build_faulted_bus_matrix(Y_pre_phys, FAULT_BUS_ID)
        Yint_fault = kron_reduce_internal(Y_phys_fault, gen_buses, xd_prime, load_adm, 
                                          is_bus_fault_active=True, fault_bus_idx=FAULT_BUS_ID)
        
    elif FAULT_TYPE == 'LINE':
        # Split Line and Ground Middle (Matrix stays 9x9)
        Y_phys_fault = build_faulted_line_matrix(Y_pre_phys, FAULT_LINE_FROM, FAULT_LINE_TO, y_line_dyn, FAULT_RATIO)
        Yint_fault = kron_reduce_internal(Y_phys_fault, gen_buses, xd_prime, load_adm, is_bus_fault_active=False)

    # ==========================================================================
    # STEP 3: POST-FAULT MATRIX (Y_post)
    # ==========================================================================
    print("Calculating Yint_post...")

    if CLEARING_ACTION == 'TRIP':
        # Line is removed from the system
        Y_phys_post = build_tripped_line_matrix(Y_pre_phys, FAULT_LINE_FROM, FAULT_LINE_TO, y_line_dyn, b_line_dyn)
        Yint_post = kron_reduce_internal(Y_phys_post, gen_buses, xd_prime, load_adm, is_bus_fault_active=False)
        
        # Save physical post-fault for Load Flow checking
        Y_phys_post_l = add_loads_to_Ybus(Y_phys_post, load_adm)

    elif CLEARING_ACTION == 'NONE':
        # System returns to Pre-Fault state
        Yint_post = Yint_pre.copy()
        Y_phys_post = Y_pre_phys.copy() 

    # ==========================================================================
    # VERIFICATION PRINT
    # ==========================================================================
    print("\n--- Python Result: Reduced Y Matrices ---")
    
    # We use slicing [:3, :3] to show only the top-left 3x3 corner for readability
    print("Yint_pre (Top 3x3):")
    print(Yint_pre[:3, :3])
    
    print("\nYint_fault (Top 3x3):")
    print(Yint_fault[:3, :3])
    
    print("\nYint_post (Top 3x3):")
    print(Yint_post[:3, :3])

    # ==========================================================================
    # PREPARE FINAL MATRICES (WITH LOADS ADDED)
    # ==========================================================================
    # 1. Create Base Loaded Matrix
    Y_pre_loaded = add_loads_to_Ybus(Y_pre_phys, load_adm)

    # 2. Create Faulted Matrix (With Loads)
    if FAULT_TYPE == 'BUS':
        # Logic: Take loaded matrix -> Remove row/col
        Y_fault_loaded = build_faulted_bus_matrix(Y_pre_loaded, FAULT_BUS_ID)
    else: # LINE
        # Logic: Take loaded matrix -> Split line
        Y_fault_loaded = build_faulted_line_matrix(Y_pre_loaded, FAULT_LINE_FROM, FAULT_LINE_TO, y_line_dyn, FAULT_RATIO)

    # 3. Create Post-Fault Matrix (With Loads)
    if CLEARING_ACTION == 'TRIP':
        # Logic: Take loaded matrix -> Trip line
        Y_post_loaded = build_tripped_line_matrix(Y_pre_loaded, FAULT_LINE_FROM, FAULT_LINE_TO, y_line_dyn, b_line_dyn)
    else: # NONE
        Y_post_loaded = Y_pre_loaded.copy()

    # ==========================================================================
    # SAVE RESULTS
    # ==========================================================================
    print(f"\nSaving to {mat_file} ...")
    
    mdict = {
        "Yi_pre":     Y_pre_loaded,
        "Y_fault":    Y_fault_loaded,
        "Y_post":     Y_post_loaded,
        "Yint_pre":   Yint_pre,
        "Yint_fault": Yint_fault,
        "Yint_post":  Yint_post
    }
    
    sio.savemat(mat_file, mdict)
    print("SUCCESS: Matrices saved.")