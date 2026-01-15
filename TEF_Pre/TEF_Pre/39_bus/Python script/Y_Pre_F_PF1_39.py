import os
import re
import numpy as np
import scipy.io as sio
import sys
import config

# Attempt to import PSS/E modules
try:
    import psse3603 as psse_mod
    import psspy
    import dyntools
except ImportError:
    pass 

# ----------------- USER CONFIG (centralized) -----------------
work_dir = config.WORK_DIR
result_dir = config.RESULT_DIR
raw_file = config.raw_path("IEEE39bus1.raw")
dyr_file = config.DYR_FILE
txt_file = config.out_path("Ybus_Export.txt")
mat_file = config.out_path("Y_all.mat")
mat_info = os.path.join(config.RESULT_DIR, "Fault_Info.mat")

# Fault & clearing specification (IEEE 39 Bus)
fault_bus = config.FAULT_BUS
trip_line_from = config.TRIP_LINE_FROM
trip_line_to = config.TRIP_LINE_TO
circuit_id = config.CKT_ID

gen_buses = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39]  # Generator bus numbers
xd_prime = [0.025, 0.05, 0.045, 0.035, 0.089, 0.04, 0.044, 0.045, 0.045, 0.004]
load_adm = {
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
# ------------------------------------------------

def init_psse_and_read_case():
    psspy.psseinit(config.PSSE_INIT)
    ierr = psspy.read(0, raw_file)
    if ierr != 0:
        raise RuntimeError(f"psspy.read returned error {ierr}")

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

def build_faulted_physical(Y_full, fault_bus_num):
    """ Removes row/col of the fault bus """
    idx = fault_bus_num - 1
    Yf = np.delete(Y_full, idx, axis=0)
    Yf = np.delete(Yf, idx, axis=1)
    return Yf

def build_postfault_physical_exact(Y_full, a_bus, b_bus, y_line, b_charging):
    """ Removes line a-b using exact extracted parameters """
    idx_a = a_bus - 1
    idx_b = b_bus - 1
    Yp = Y_full.copy()
    
    # Subtract Y_series and B_charging from diagonals
    y_total = y_line + b_charging
    Yp[idx_a, idx_a] -= y_total
    Yp[idx_b, idx_b] -= y_total
    
    # Open the line (Zero off-diagonals)
    Yp[idx_a, idx_b] = 0.0 + 0.0j
    Yp[idx_b, idx_a] = 0.0 + 0.0j
    return Yp

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
    
    # -------------------------------------------------------------
    # FIX: Check against 39 (Total Buses) instead of 9
    # If matrix size N is 38 (less than 39), it means a bus (fault bus) was removed.
    # -------------------------------------------------------------
    is_faulted = (N < 39) 

    # Add Gen Admittances to Diagonals
    gen_indices_map = [] # Stores (k, physical_index)
    for i, bus in enumerate(gen_buses_list):
        if is_faulted:
            if bus == fault_bus: continue 
            # Index shift: if gen bus is after fault bus, subtract 1 (bus-2 for 0-based index)
            phys_idx = bus - 1 if bus < fault_bus else bus - 2
        else:
            # Normal indexing (bus-1 for 0-based index)
            phys_idx = bus - 1
            
        gen_indices_map.append((i, phys_idx))
        YD[phys_idx, phys_idx] += ybar[i, i]

    # Add Load Admittances to Diagonals
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

# === Main execution ===
if __name__ == "__main__":
    if not os.path.exists(raw_file):
        raise FileNotFoundError(f"raw file not found: {raw_file}")
    
    print("Initializing PSS/E...")
    init_psse_and_read_case()

    # --- NEW: Extract Line Parameters Automatically ---
    print(f"Extracting parameters for line {trip_line_from}-{trip_line_to}...")
    Y_line_dynamic, B_line_dynamic = get_line_params_from_psse(trip_line_from, trip_line_to, circuit_id)
    # --------------------------------------------------

    print("Solving (No Load Conversion)...")
    solve_and_convert() 

    print("Exporting Ybus...")
    export_ybus_to_text()
    
    print("Parsing Ybus...")
    Y_pre = parse_ybus_text_to_numpy()
    
    # We typically don't need Y_pre_l for the function calls below 
    # because 'kron_reduce_internal' adds the loads internally via `load_admittance`.
    # However, kept here if you need to debug the full loaded matrix.
    Y_pre_l = add_loads_to_Ybus(Y_pre, load_adm)
    
    # -------------------------------------------------------------
    # 1. PRE-FAULT (With Loads Added Manually)
    # -------------------------------------------------------------
    print("Calculating Pre-Fault Reduced Matrix...")
    Yint_pre = kron_reduce_internal(Y_pre, gen_buses, xd_prime, load_admittance=load_adm)

    # -------------------------------------------------------------
    # 2. FAULTED (With Loads Added Manually)
    # -------------------------------------------------------------
    print(f"Calculating Faulted Reduced Matrix...")
    Y_fault = build_faulted_physical(Y_pre, fault_bus)
    Yint_fault = kron_reduce_internal(Y_fault, gen_buses, xd_prime, load_admittance=load_adm)

    # -------------------------------------------------------------
    # 3. POST-FAULT (With Loads Added Manually)
    # -------------------------------------------------------------
    print(f"Calculating Post-Fault Reduced Matrix...")
    
    # Using the automatically extracted Y_line and B_line here:
    Y_post = build_postfault_physical_exact(
        Y_pre, 
        trip_line_to, 
        trip_line_from, 
        Y_line_dynamic, 
        B_line_dynamic
    )
    Y_post1= build_postfault_physical_exact(
        Y_pre_l, 
        trip_line_from, 
        trip_line_to, 
        Y_line_dynamic, 
        B_line_dynamic
    )
    
    Yint_post = kron_reduce_internal(Y_post, gen_buses, xd_prime, load_admittance=load_adm)

    # Verification Print
    print("\n--- Python Result: Reduced Y Matrices ---")
    print("Yint_pre (Top 10x10):")
    print(Yint_pre)
    print("\nYint_fault (Top 10x10):")
    print(Yint_fault)
    print("\nYint_post (Top 10x10):")
    print(Yint_post)
    
    # Save
    print(f"\nSaving to {mat_file} ...")
    mdict = {
        "Yi_pre": Y_pre,
        "Y_fault": Y_fault,
        "Y_post": Y_pre,
        "Yint_pre": Yint_pre,
        "Yint_fault": Yint_fault,
        "Yint_post": Yint_pre
    }
    sio.savemat(mat_file, mdict)
     
    sio.savemat(mat_info, {
        "Fault_Type": "Line_Fault",
        "Fault_Location": fault_bus,
        "Tripped_Line": [trip_line_from, trip_line_to],
        "Fault_Impedance": 0
    })
    print("SUCCESS: Matrices saved.")

    