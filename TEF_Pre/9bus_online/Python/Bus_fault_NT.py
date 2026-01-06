# script_1_bus_fault.py
import os
import re
import numpy as np
import scipy.io as sio
import psse3603 # type: ignore
import psspy    # type: ignore
import config

# --- SETUP ---
raw_file = config.RAW_FILE
txt_file = config.TXT_FILE_YBUS
mat_file = config.MAT_FILE_YBUS
mat_info = os.path.join(config.RESULT_DIR, "Fault_Info.mat")

fault_bus = config.FAULT_BUS
gen_buses = config.GEN_BUSES
xd_prime = config.XD_PRIME
load_adm = config.LOAD_ADM

# --- HELPER FUNCTIONS ---
def init_psse():
    psspy.psseinit(200)
    psspy.read(0, raw_file)

def solve_and_convert():
    psspy.fnsl([0,0,0,1,1,0,99,0])
    psspy.cong(0)
    try:
        psspy.ordr(0)
        psspy.fact()
        psspy.tysl(0)
    except: pass

def export_and_parse_ybus():
    psspy.output_y_matrix(0, 1, 0, 0, txt_file)
    Yreal, Yimag = {}, {}
    max_bus = 0
    float_re = re.compile(r"[-+]?\d*\.\d+|\d+")
    with open(txt_file, "r") as f:
        for line in f:
            parts = float_re.findall(line)
            if len(parts) >= 4:
                try:
                    i, j = int(parts[0]), int(parts[1])
                    val = float(parts[2]) + 1j*float(parts[3])
                    max_bus = max(max_bus, i, j)
                    Yreal[(i, j)] = val.real
                    Yimag[(i, j)] = val.imag
                except: continue
    Y = np.zeros((max_bus, max_bus), dtype=complex)
    for (i, j), r in Yreal.items():
        im = Yimag.get((i, j), 0.0)
        Y[i-1, j-1] = r + 1j*im
        Y[j-1, i-1] = r + 1j*im
    return Y

def add_loads(Y, loads):
    Y_new = Y.copy()
    for b, yl in loads.items():
        Y_new[b-1, b-1] += yl
    return Y_new

def kron_reduce(Y_full, remove_bus=None):
    # If Y_full already has loads, we don't need to add them again if we passed the loaded matrix
    # But for safety in this function, let's assume Y_full IS the system matrix to be reduced.
    # We just handle generator nodes.
    
    Y_sys = Y_full.copy()
    
    # Handle Fault Removal
    if remove_bus:
        idx = remove_bus - 1
        Y_sys = np.delete(Y_sys, idx, 0)
        Y_sys = np.delete(Y_sys, idx, 1)

    K = len(gen_buses)
    N = Y_sys.shape[0]
    
    ybar = np.zeros((K, K), dtype=complex)
    YD = Y_sys.copy()
    gen_map = []
    
    for i, gb in enumerate(gen_buses):
        ybar[i,i] = 1.0/(1j * xd_prime[i])
        if remove_bus:
            pidx = gb - 1 if gb < remove_bus else gb - 2
        else:
            pidx = gb - 1
        gen_map.append((i, pidx))
        YD[pidx, pidx] += ybar[i,i]
        
    YA = ybar.copy()
    YB = np.zeros((K, N), dtype=complex)
    for k, pidx in gen_map:
        YB[k, pidx] = -ybar[k,k]
    YC = YB.T
    
    try:
        Y_red = YA - np.dot(YB, np.linalg.solve(YD, YC))
    except:
        Y_red = np.zeros((K, K))
    return Y_red

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    init_psse()
    solve_and_convert()
    Y_pre = export_and_parse_ybus()          # Raw PSS/E (Lines/Shunts/Xfmrs)
    Y_pre_l = add_loads(Y_pre, load_adm)     # Pre-Fault Physical (With Loads)
    
    print("--- GENERATING BUS FAULT (NO TRIP) ---")
    
    # 1. Pre-Fault Reduced
    # Note: We pass Y_pre_l so loads are included in reduction
    Yint_pre = kron_reduce(Y_pre_l)
    
    # 2. Faulted
    # Physical: Remove row/col of fault bus from Y_pre (the raw matrix usually, or loaded? 
    # Usually fault current dominates load, but strictly: Y_pre_l with row/col removed)
    Y_fault_phys = np.delete(Y_pre, fault_bus-1, 0)
    Y_fault_phys = np.delete(Y_fault_phys, fault_bus-1, 1)
    
    # Reduced:
    Yint_fault = kron_reduce(Y_pre_l, remove_bus=fault_bus)
    
    # 3. Post-Fault (No Trip)
    # Physical: Same as Pre-Fault Loaded
    Y_post1 = Y_pre_l.copy()
    
    # Reduced: Same as Pre-Fault Reduced
    Yint_post = kron_reduce(Y_post1)
    
    # SAVE
    mdict = {
        "Yi_pre": Y_pre,           # Raw PSS/E
        "Y_fault": Y_fault_phys,   # Unreduced Fault
        "Y_post": Y_post1,         # Unreduced Post-Fault (Loaded)
        "Yint_pre": Yint_pre,
        "Yint_fault": Yint_fault,
        "Yint_post": Yint_post
    }
    sio.savemat(mat_file, mdict)
    
    # Info
    sio.savemat(mat_info, {
        "Fault_Type": "Bus_Fault",
        "Fault_Location": fault_bus,
        "Tripped_Line": [],
        "Fault_Impedance": 0
    })
    print(f"Saved Bus Fault data to {mat_file}")