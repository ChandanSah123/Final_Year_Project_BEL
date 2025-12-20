import os
import re
import numpy as np
import scipy.io as sio

# Attempt to import PSS/E modules
try:
    import psse3603 as psse_mod
    import psspy
    import dyntools
except ImportError:
    pass # Will fail later if not available

# ----------------- USER CONFIG -----------------
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_9 bus\Python Script"
result_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_9 bus\Matlab Script"
raw_file = os.path.join(work_dir, "IEEE9bus.raw")
dyr_file = os.path.join(work_dir, "ieee9bus.dyr")
txt_file = os.path.join(result_dir, "Ybus_Export.txt")
mat_file = os.path.join(result_dir, "Y_all.mat")

# Fault & clearing specification
fault_bus = 7
trip_line_from = 7
trip_line_to = 5

gen_buses = [1, 2, 3]
xd_prime = [0.0608, 0.1198, 0.1813]

# Hardcoded Load Admittances (Exact match to MATLAB)
load_adm = {
    5: 1.26 - 0.504j,
    6: 0.877 - 0.292j,
    8: 0.969 - 0.339j
}

# Exact Post-Fault Line Parameters (From MATLAB)
# Used to remove the line 5-7
Y_line_57 = 1.187 - 5.975j
B_line_57 = 0.153j

# ------------------------------------------------

def init_psse_and_read_case():
    psspy.psseinit(200)
    ierr = psspy.read(0, raw_file)
    if ierr != 0:
        raise RuntimeError(f"psspy.read returned error {ierr}")

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
    """ Removes line a-b using exact MATLAB parameters """
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
    
    # Detect if we are in Faulted Mode (Bus removed)
    # IEEE 9 bus is size 9. If size is 8, it's faulted.
    is_faulted = (N < 9) 

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
def add_loads_to_Ybus(Ybus,load_admittance):
    Y_loaded=Ybus.copy()
    for bus_num,yload in load_admittance.items():
        idx=bus_num-1
        Y_loaded[idx,idx]+=yload
        return Y_loaded

# === Main execution ===
if __name__ == "__main__":
    if not os.path.exists(raw_file):
        raise FileNotFoundError(f"raw file not found: {raw_file}")
    
    print("Initializing PSS/E...")
    init_psse_and_read_case()

    print("Solving (No Load Conversion)...")
    solve_and_convert() 

    print("Exporting Ybus...")
    export_ybus_to_text()
    
    print("Parsing Ybus...")
    Y_pre = parse_ybus_text_to_numpy()
    Y_pre_l=add_loads_to_Ybus(Y_pre,load_adm)
    
    # -------------------------------------------------------------
    # 1. PRE-FAULT (With Loads Added Manually)
    # -------------------------------------------------------------
    print("Calculating Pre-Fault Reduced Matrix...")
    # PASS load_adm HERE!
    Yint_pre = kron_reduce_internal(Y_pre, gen_buses, xd_prime, load_admittance=load_adm)

    # -------------------------------------------------------------
    # 2. FAULTED (With Loads Added Manually)
    # -------------------------------------------------------------
    print(f"Calculating Faulted Reduced Matrix...")
    Y_fault = build_faulted_physical(Y_pre, fault_bus)
    Y_fault_l = build_faulted_physical(Y_pre_l, fault_bus)
    # PASS load_adm HERE!
    Yint_fault = kron_reduce_internal(Y_fault, gen_buses, xd_prime, load_admittance=load_adm)

    # -------------------------------------------------------------
    # 3. POST-FAULT (With Loads Added Manually)
    # -------------------------------------------------------------
    print(f"Calculating Post-Fault Reduced Matrix...")
    Y_post = build_postfault_physical_exact(Y_pre, trip_line_to, trip_line_from, Y_line_57, B_line_57)
    Y_post_l=add_loads_to_Ybus(Y_post,load_adm)
    # PASS load_adm HERE!
    Yint_post = kron_reduce_internal(Y_post, gen_buses, xd_prime, load_admittance=load_adm)

    # Verification Print
    print("\n--- Python Result: Yint_Fault (Top 3x3) ---")
    print("\n--- Python Result:  ---")
    print("Yint_pre")
    print(Yint_pre)
    print("Yint_fault")
    print(Yint_fault)
    print("Yint_post")
    print(Yint_post)
    
    # Save
    print(f"\nSaving to {mat_file} ...")
    mdict = {
      
        "Yint_pre": Yint_pre,
        "Yint_fault": Yint_fault,
        "Yint_post": Yint_post
    }
    sio.savemat(mat_file, mdict)
    print("SUCCESS: Matrices saved.")