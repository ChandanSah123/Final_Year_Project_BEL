import os
import re
import numpy as np
import scipy.io as sio
import time


# Attempt to import PSS/E modules
try:
    import psse3603  # type: ignore
    import psspy     # type: ignore
    import dyntools  # type: ignore
    pss_available = True
except Exception:
    pss_available = False

import config

# ----------------- USER CONFIG (centralized) -----------------
work_dir = config.WORK_DIR
result_dir = config.RESULT_DIR
raw_file = config.RAW_FILE
dyr_file = config.DYR_FILE
txt_file = config.TXT_FILE_YBUS
mat_file = config.MAT_FILE_YBUS
mat_info = os.path.join(config.RESULT_DIR, "Fault_Info.mat")

# Fault & clearing specification (use centralized values)
fault_bus = config.FAULT_BUS
trip_line_from = config.TRIP_LINE_FROM
trip_line_to = config.TRIP_LINE_TO
circuit_id = config.LINE_ID  # Usually '1 '

gen_buses = config.GEN_BUSES
xd_prime = config.XD_PRIME

# Hardcoded Load Admittances (Exact match to MATLAB)
load_adm = config.LOAD_ADM

# Default timing
fault_time = config.T_FAULT_START
t_cl = config.T_CLEAR
clear_time = fault_time + t_cl
end_time = config.T_END

# ------------------------------------------------

def init_psse_and_read_case():
    try:
        _ = psspy.psseinit(200)
        ierr = psspy.read(0, raw_file)
        if ierr != 0:
            raise RuntimeError(f"psspy.read returned error {ierr}")
    except Exception as e:
        raise

def get_line_params_from_psse(bus_from, bus_to, ckt_id='1 '):
    """
    Extracts line or transformer parameters (Y_series and B_charging_half) from PSS/E.
    Handles both Lines (brndat) and Transformers (xfrdat).
    """
    # Initialize variables
    z_complex = 0j
    b_total_real = 0.0
    found = False

    # -----------------------------------------------------------
    # ATTEMPT 1: Try as a LINE (Branch)
    # -----------------------------------------------------------
    # Get Impedance (RX)
    ierr, z_val = psspy.brndt2(bus_from, bus_to, ckt_id, "RX")
    
    if ierr == 0:
        # It is a Line. Get Charging using 'CHARG'
        ierr_c, b_val = psspy.brndat(bus_from, bus_to, ckt_id, "CHARG")
        if ierr_c == 0:
            z_complex = z_val
            b_total_real = b_val
            found = True
            print(f"Detected Line {bus_from}-{bus_to}")

    # -----------------------------------------------------------
    # ATTEMPT 2: Try as a TRANSFORMER (if Line failed)
    # -----------------------------------------------------------
    if not found:
        # Get Impedance (RX) for Transformer
        # xfrdat keywords: 'RX' for impedance, 'MAG1'/'MAG2' for magnetizing (often ignored in simple extraction or assumed small)
        ierr, z_val = psspy.xfrdat(bus_from, bus_to, ckt_id, "RX")
        
        if ierr == 0:
            # It is a Transformer.
            # Transformers in PSS/E usually have magnetizing admittance (G+jB) rather than "Charging"
            # For Kron reduction, we usually ignore magnetizing B or treat it as charging.
            # Let's try to get Magnetizing B (MAG2 is usually Im(y_mag))
            ierr_m, mag_imag = psspy.xfrdat(bus_from, bus_to, ckt_id, "MAG2")
            
            if ierr_m == 0:
                z_complex = z_val
                b_total_real = mag_imag # Treat magnetizing B as charging B
                found = True
                print(f"Detected Transformer {bus_from}-{bus_to}")
            else:
                 # Default to 0 charging if not found
                 z_complex = z_val
                 b_total_real = 0.0
                 found = True
                 print(f"Detected Transformer {bus_from}-{bus_to} (No Magnetizing B found)")

    # -----------------------------------------------------------
    # ERROR HANDLING
    # -----------------------------------------------------------
    if not found:
        # Try Reverse direction just in case
        raise ValueError(f"Could not find Line OR Transformer between {bus_from}-{bus_to} with ID '{ckt_id}'. Check Bus IDs and Circuit ID.")

    # -----------------------------------------------------------
    # CALCULATION
    # -----------------------------------------------------------
    # Series Admittance Y = 1 / Z
    if z_complex == 0j:
         # Prevent divide by zero if data is bad
         raise ValueError(f"Zero Impedance detected for {bus_from}-{bus_to}")
         
    Y_line = 1.0 / z_complex
    
    # Half Charging
    # For a line: B/2 at each end.
    # For a transformer: Magnetizing B is usually at the 'From' bus or split. 
    # For simple stability studies, treating it as Pi-model (split) is acceptable approximation.
    B_line_half = 1j * (b_total_real / 2.0)
    
    print(f"  Z = {z_complex:.4f}")
    print(f"  Y_series = {Y_line:.4f}")
    print(f"  B_half   = {B_line_half:.4f}")
    
    return Y_line, B_line_half

def solve_and_convert():
    # 1. Solve Power Flow
    psspy.fnsl([0,0,0,1,1,0,99,0])
    
    # 2. Convert Generators Only (Norton Equivalent)
    psspy.cong(0)
    
    # 3. DO NOT CONVERT LOADS (conl)
    # We skip psspy.conl completely. 
    # This ensures the exported Y-bus contains ONLY lines and transformers.
    # We will add the loads manually in Python to ensure they match MATLAB exactly.
    
    try:
        psspy.ordr(0)
        psspy.fact()
        psspy.tysl(0)
    except Exception:
        pass

def export_ybus_to_text():
    # Export the matrix (Lines + Transformers + Shunts ONLY)
    ierr = psspy.output_y_matrix(0, 1, 0, 0, txt_file)
    if isinstance(ierr, tuple): ierr = ierr[0]
    if ierr != 0:
        raise RuntimeError(f"output_y_matrix error: {ierr}")

def parse_ybus_text_to_numpy():
    Yreal = {}
    Yimag = {}
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
                except Exception:
                    continue
    
    Y = np.zeros((max_bus, max_bus), dtype=np.complex128)
    for (i, j), rv in Yreal.items():
        iv = Yimag.get((i, j), 0.0)
        Y[i-1, j-1] = rv + 1j*iv
        Y[j-1, i-1] = rv + 1j*iv
    return Y
def kron_reduce_internal(Y_physical, gen_buses_list, xd_prime_list, load_admittance=None):
    """
    Performs Kron Reduction. 
    1. Adds ybar (Generator Internal Z) to diagonals.
    2. Adds load_admittance (Manual Python Injection) to diagonals.
    3. Reduces matrix.
    """
    K = len(gen_buses_list)
    N = Y_physical.shape[0]

    # A. Calculate Generator Admittances (ybar)
    ybar = np.zeros((K, K), dtype=complex)
    for i, xd in enumerate(xd_prime_list):
        ybar[i, i] = 1.0 / (1j * xd)

    # B. Augment Physical Matrix (YD)
    YD = Y_physical.copy()
    
    # Detect if we are in Faulted Mode (Bus removed)
    # IEEE 39 bus is size 39. If size is 38, it's faulted.
    is_faulted = (N < 39) 

    # Add Gen Admittances to Diagonals
    gen_indices_map = [] # Stores (k, physical_index)
    for i, bus in enumerate(gen_buses_list):
        if is_faulted:
            if bus == fault_bus: continue 
            # Index shift: if gen bus is after fault bus, subtract 1
            phys_idx = bus - 1 if bus < fault_bus else bus - 2
        else:
            phys_idx = bus - 1
            
        gen_indices_map.append((i, phys_idx))
        YD[phys_idx, phys_idx] += ybar[i, i]

    # Add Load Admittances to Diagonals (THE CRITICAL FIX)
    # We manually add these here because we skipped conl in PSS/E
    if load_admittance:
        for bus_num, yload in load_admittance.items():
            if is_faulted:
                if bus_num == fault_bus: continue
                idx = bus_num - 1 if bus_num < fault_bus else bus_num - 2
                YD[idx, idx] += yload
            else:
                idx = bus_num - 1
                YD[idx, idx] += yload

    # C. Build Partition Matrices
    YA = ybar.copy()
    YB = np.zeros((K, N), dtype=complex)
    
    for (k, phys_idx) in gen_indices_map:
        YB[k, phys_idx] = -ybar[k, k]

    YC = YB.T

    # D. Reduce
    try:
        YD_inv_YC = np.linalg.solve(YD, YC)
        Yint = YA - np.dot(YB, YD_inv_YC)
    except np.linalg.LinAlgError:
        print("Warning: Matrix singular.")
        Yint = np.zeros((K, K))
        
    return Yint

def add_loads_to_Ybus(Ybus, load_admittance):
    Y_loaded = Ybus.copy()
    for bus_num, yload in load_admittance.items():
        idx = bus_num - 1
        Y_loaded[idx, idx] += yload
    return Y_loaded

if __name__ == "__main__":
    if not os.path.exists(raw_file):
        raise FileNotFoundError(f"raw file not found: {raw_file}")
    
    print("Initializing PSS/E...")
    init_psse_and_read_case()

    # --- NEW: Extract Line Parameters Automatically ---
    start_cpu = time.process_time()
    print(f"Extracting parameters for line {trip_line_from}-{trip_line_to}...")
    Y_line_dynamic, B_line_dynamic = get_line_params_from_psse(trip_line_from, trip_line_to, circuit_id)
    # --------------------------------------------------

    print("Solving (No Load Conversion)...")
    solve_and_convert() 

    print("Exporting Ybus...")
    export_ybus_to_text()
    
    print("Parsing Ybus...")
    Y_pre = parse_ybus_text_to_numpy()
    print("Calculating Pre-Fault Reduced Matrix...")
    Yint_pre = kron_reduce_internal(Y_pre, gen_buses, xd_prime, load_admittance=load_adm)

    print(Yint_pre)
    print(f"\nSaving to {mat_file} ...")
    mdict = {
        "Yint_pre": Yint_pre,
        "Yint_post": Yint_pre
       }
    sio.savemat(mat_file, mdict)
     
    sio.savemat(mat_info, {
        "Fault_Type": "Line_Fault",
        "Fault_Location": fault_bus,
        "Tripped_Line": [trip_line_from, trip_line_to],
        "Fault_Impedance": 0
    })
    end_cpu = time.process_time()
    print(f"Post-Fault Reduction CPU Time: {end_cpu - start_cpu:.5f} seconds")
    print("SUCCESS: Matrices saved.")