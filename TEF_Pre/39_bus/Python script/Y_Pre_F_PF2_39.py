import os
import re
import numpy as np
import scipy.io as sio
import sys

# Attempt to import PSS/E modules
try:
    import psse3603 as psse_mod
    import psspy
    import dyntools
except ImportError:
    pass 

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================
# --- File Paths ---
work_dir = r"C:\Users\Acer\Desktop\final year project\Final_Year_Project_BEL\TEF_Pre\39_bus\Python script"
result_dir = r"C:\Users\Acer\Desktop\final year project\Final_Year_Project_BEL\TEF_Pre\39_bus\Matlab script"

raw_file   = os.path.join(work_dir, "IEEE39bus.raw")    # Ensure this matches your filename
dyr_file   = os.path.join(work_dir, "ieee39buscls.dyr") # Your GENCLS dyr file
txt_file   = os.path.join(result_dir, "Ybus_Export.txt")
mat_file   = os.path.join(result_dir, "Y_all.mat")

# --- Scenario Settings ---
# Fault Specification
FAULT_BUS      = 16       # The bus where the 3-phase fault occurs
TRIP_LINE_FROM = 16       # Start of the line to trip
TRIP_LINE_TO   = 17       # End of the line to trip
CKT_ID         = '1 '     # Circuit ID of the line to trip

# ==============================================================================
# 2. DATA EXTRACTION FUNCTIONS
# ==============================================================================

def init_psse_and_read_case():
    psspy.psseinit(2000) # Support larger systems
    
    # Read RAW
    if not os.path.exists(raw_file):
        raise FileNotFoundError(f"RAW file not found at: {raw_file}")
    ierr = psspy.read(0, raw_file)
    if ierr != 0: raise RuntimeError(f"Error reading RAW file: {ierr}")
    
    # Read DYR (To apply the GENCLS models to the buses)
    if os.path.exists(dyr_file):
        print(f"Reading DYR file: {dyr_file}")
        psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")
    else:
        print("Warning: DYR file not found. Assuming GENCLS models are implied or handled manually.")

def extract_machine_data_gencls():
    """
    Specific extraction for GENCLS.
    GENCLS uses the Generator Source Impedance (ZSORCE) from the RAW file 
    as the internal transient reactance (X'd).
    """
    print("Extracting GENCLS parameters (Z_source)...")
    
    # Get all machine buses
    ierr, mach_buses = psspy.amachint(-1, 4, "NUMBER")
    if ierr != 0: raise RuntimeError("Error retrieving machines.")
    mach_buses = mach_buses[0]
    
    buses_out = []
    xd_out = []
    
    for bus in mach_buses:
        # Check if machine is in-service
        ierr, status = psspy.macint(bus, '1 ', "STATUS")
        if status == 0: continue

        # For GENCLS, we need ZSORCE (specifically the imaginary part X)
        ierr, z_val = psspy.macdat(bus, '1 ', "ZSORCE")
        
        if ierr == 0:
            # z_val is a complex number (R + jX)
            # We strictly need X for the reactance matrix
            x_val = z_val.imag
            
            # Sanity check: If X is 0, it might be an error or infinite bus
            if x_val < 0.0001:
                print(f"  Warning: Bus {bus} has X_source ~ 0. Using small epsilon 0.0001")
                x_val = 0.0001
                
            buses_out.append(bus)
            xd_out.append(x_val)
            print(f"  Bus {bus}: Found GENCLS/Machine. X'd (Zsource) = {x_val:.5f}")
        else:
            print(f"  Error reading ZSORCE for Bus {bus}")

    # Sort by bus number for consistency
    sorted_pairs = sorted(zip(buses_out, xd_out))
    return [p[0] for p in sorted_pairs], [p[1] for p in sorted_pairs]

def extract_load_admittances():
    """ Calculates Y_load = (P - jQ) / |V|^2 from the solved state. """
    load_map = {}
    ierr, buses = psspy.abusint(-1, 2, "NUMBER")
    buses = buses[0]
    
    for bus in buses:
        # Get Total Load at Bus
        ierr, cmpval = psspy.loddt2(bus, '1 ', "MVA", "ACT")
        if ierr == 0:
            p_load = cmpval.real
            q_load = cmpval.imag
            
            if abs(p_load) > 1e-6 or abs(q_load) > 1e-6:
                ierr_v, v_mag = psspy.busdat(bus, "PU")
                if ierr_v == 0 and v_mag > 0.001:
                    # Y = S* / |V|^2  (S in pu, PSS/E load in MVA, Base=100)
                    s_pu = (p_load - 1j * q_load) / 100.0
                    y_load = s_pu / (v_mag ** 2)
                    load_map[bus] = y_load
    
    print(f"Extracted {len(load_map)} loads.")
    return load_map

def get_line_params(bus_from, bus_to, ckt_id='1 '):
    """ Extracts Y_series and B_charging/2 for a specific line/trafo. """
    # 1. Try Line
    ierr, z_val = psspy.brndt2(bus_from, bus_to, ckt_id, "RX")
    if ierr == 0:
        ierr_c, b_val = psspy.brndat(bus_from, bus_to, ckt_id, "CHARG")
        return 1.0/z_val, 1j*(b_val/2.0)
    
    # 2. Try Transformer
    ierr, z_val = psspy.xfrdat(bus_from, bus_to, ckt_id, "RX")
    if ierr == 0:
        ierr_m, mag_imag = psspy.xfrdat(bus_from, bus_to, ckt_id, "MAG2") 
        if ierr_m != 0: mag_imag = 0.0
        return 1.0/z_val, 1j*(mag_imag/2.0)

    raise ValueError(f"Branch {bus_from}-{bus_to} not found.")

# ==============================================================================
# 3. MATRIX MATH FUNCTIONS
# ==============================================================================

def solve_and_export_ybus():
    # Solve PF
    psspy.fnsl([0,0,0,1,1,0,99,0])
    # Generator Norton Equivalent
    psspy.cong(0) 
    # NO Load conversion (we handle loads in Python)
    try:
        psspy.ordr(0)
        psspy.fact()
        psspy.tysl(0)
    except: pass
    
    ierr = psspy.output_y_matrix(0, 1, 0, 0, txt_file)
    if ierr != 0 and ierr != (0,0): 
        print(f"Warning: Ybus export returned {ierr}")

def parse_ybus_file():
    Yreal, Yimag = {}, {}
    max_bus = 0
    float_re = re.compile(r"[-+]?\d*\.\d+|\d+")
    
    with open(txt_file, "r") as f:
        for line in f:
            parts = float_re.findall(line)
            if len(parts) >= 4:
                i, j = int(parts[0]), int(parts[1])
                real, imag = float(parts[2]), float(parts[3])
                max_bus = max(max_bus, i, j)
                Yreal[(i, j)] = real
                Yimag[(i, j)] = imag
                
    dim = max(max_bus, 39)
    Y = np.zeros((dim, dim), dtype=np.complex128)
    for (i, j), rv in Yreal.items():
        iv = Yimag.get((i, j), 0.0)
        Y[i-1, j-1] = rv + 1j*iv
        Y[j-1, i-1] = rv + 1j*iv
    return Y

def kron_reduce(Y_full, gen_buses, xd_primes, load_adm, fault_bus=None):
    N = Y_full.shape[0]
    K = len(gen_buses)
    is_fault_matrix = (fault_bus is not None)
    
    YD = Y_full.copy()
    
    def get_idx(bus_id):
        if is_fault_matrix:
            if bus_id == fault_bus: return None
            return bus_id - 1 if bus_id < fault_bus else bus_id - 2
        return bus_id - 1

    # Add Loads
    if load_adm:
        for bus, y in load_adm.items():
            idx = get_idx(bus)
            if idx is not None:
                YD[idx, idx] += y

    # Add Generators (1/jXd')
    ybar = np.zeros((K, K), dtype=complex)
    gen_map = [] 

    for k, (bus, xd) in enumerate(zip(gen_buses, xd_primes)):
        y_gen = 1.0 / (1j * xd)
        ybar[k, k] = y_gen
        
        phys_idx = get_idx(bus)
        if phys_idx is not None:
            YD[phys_idx, phys_idx] += y_gen
            gen_map.append((k, phys_idx))

    # Partitioning & Reduction
    YA = ybar 
    YB = np.zeros((K, YD.shape[0]), dtype=complex)
    for (k, phys_idx) in gen_map:
        YB[k, phys_idx] = -ybar[k, k]
    YC = YB.T
    
    try:
        Y_int = YA - np.dot(YB, np.linalg.solve(YD, YC))
    except np.linalg.LinAlgError:
        print("Error: Singular Matrix during Kron Reduction.")
        return np.zeros((K, K))
        
    return Y_int

# ==============================================================================
# 4. MAIN EXECUTION
# ==============================================================================
if __name__ == "__main__":
    print("--- GENCLS AUTOMATED YBUS EXTRACTION ---")
    
    # 1. Init & Read
    init_psse_and_read_case()
    
    # 2. Extract Data (Using GENCLS specific logic)
    gens, xds = extract_machine_data_gencls()
    print(f"  Found {len(gens)} GENCLS machines.")
    
    print("Solving PF & Extracting Loads...")
    solve_and_export_ybus()
    loads = extract_load_admittances()
    
    print(f"Extracting Line Params for {TRIP_LINE_FROM}-{TRIP_LINE_TO}...")
    y_line, b_line = get_line_params(TRIP_LINE_FROM, TRIP_LINE_TO, CKT_ID)
    
    # 3. Parse Physical Matrix
    Y_phys_base = parse_ybus_file()
    
    # --- STEP 1: PRE-FAULT ---
    print("Calculating Pre-Fault...")
    Yint_pre = kron_reduce(Y_phys_base, gens, xds, loads)
    
    # --- STEP 2: FAULTED ---
    print(f"Calculating Faulted (Bus {FAULT_BUS})...")
    idx_f = FAULT_BUS - 1
    Y_phys_fault = np.delete(Y_phys_base, idx_f, axis=0)
    Y_phys_fault = np.delete(Y_phys_fault, idx_f, axis=1)
    Yint_fault = kron_reduce(Y_phys_fault, gens, xds, loads, fault_bus=FAULT_BUS)
    
    # --- STEP 3: POST-FAULT ---
    print("Calculating Post-Fault...")
    idx_a, idx_b = TRIP_LINE_FROM - 1, TRIP_LINE_TO - 1
    Y_phys_post = Y_phys_base.copy()
    
    y_total = y_line + b_line
    Y_phys_post[idx_a, idx_a] -= y_total
    Y_phys_post[idx_b, idx_b] -= y_total
    Y_phys_post[idx_a, idx_b] = 0j
    Y_phys_post[idx_b, idx_a] = 0j
    
    Yint_post = kron_reduce(Y_phys_post, gens, xds, loads)
    
    # --- SAVE ---
    def add_loads(Y, lds):
        Yn = Y.copy()
        for b, y in lds.items():
            if b <= Yn.shape[0]: Yn[b-1, b-1] += y
        return Yn

    Y_pre_loaded = add_loads(Y_phys_base, loads)
    Y_post_loaded = add_loads(Y_phys_post, loads)

    print(f"Saving to {mat_file}...")
    sio.savemat(mat_file, {
        "Yi_pre":     Y_pre_loaded,
        "Y_post":     Y_post_loaded,
        "Yint_pre":   Yint_pre,
        "Yint_fault": Yint_fault,
        "Yint_post":  Yint_post,
        "gen_buses":  gens,
        "xd_primes":  xds
    })
    print("Done.")