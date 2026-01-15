# script_3_midline_fault.py
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

trip_from = config.TRIP_LINE_FROM
trip_to = config.TRIP_LINE_TO
ckt_id = config.LINE_ID
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
        psspy.ordr(0); psspy.fact(); psspy.tysl(0)
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

def kron_reduce(Y_full):
    # Standard Kron Reduction (No bus removal supported/needed here because 
    # we modify the matrix values directly for the fault)
    K = len(gen_buses)
    N = Y_full.shape[0]
    ybar = np.zeros((K, K), dtype=complex)
    YD = Y_full.copy()
    gen_map = []
    for i, gb in enumerate(gen_buses):
        ybar[i,i] = 1.0/(1j * xd_prime[i])
        pidx = gb - 1
        gen_map.append((i, pidx))
        YD[pidx, pidx] += ybar[i,i]
    YA = ybar.copy()
    YB = np.zeros((K, N), dtype=complex)
    for k, pidx in gen_map: YB[k, pidx] = -ybar[k,k]
    YC = YB.T
    try: Y_red = YA - np.dot(YB, np.linalg.solve(YD, YC))
    except: Y_red = np.zeros((K, K))
    return Y_red

def get_line_params(i, j, ckt):
    ierr, z = psspy.brndt2(i, j, ckt, "RX")
    ierr, chg = psspy.brndat(i, j, ckt, "CHARG")
    return 1.0/z, 1j*(chg/2.0)

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    init_psse()
    Y_line, B_half = get_line_params(trip_from, trip_to, ckt_id)
    solve_and_convert()
    Y_pre = export_and_parse_ybus()
    Y_pre_l = add_loads(Y_pre, load_adm)
    
    print("--- GENERATING 50% LINE FAULT ---")
    
    # 1. Pre-Fault Reduced
    Yint_pre = kron_reduce(Y_pre_l)
    
    # 2. Faulted (Physical Construction)
    # Start with Loaded Pre-Fault
    Y_fault_phys = Y_pre_l.copy()
    idx_f, idx_t = trip_from - 1, trip_to - 1
    
    # A. Remove the healthy line
    y_orig_total = Y_line + B_half
    Y_fault_phys[idx_f, idx_f] -= y_orig_total
    Y_fault_phys[idx_t, idx_t] -= y_orig_total
    Y_fault_phys[idx_f, idx_t] = 0.0
    Y_fault_phys[idx_t, idx_f] = 0.0
    
    # B. Add Split Shunts (Fault at 50%)
    # Z_half = Z/2 => Y_half = 2*Y
    # B_half_segment = B_total/4 (Half of the half-charging)
    # The actual B_half variable is B_total/2. So we want B_half/2.
    y_split_shunt = (2 * Y_line) + (B_half / 2.0)
    
    Y_fault_phys[idx_f, idx_f] += y_split_shunt
    Y_fault_phys[idx_t, idx_t] += y_split_shunt
    
    # Reduced Fault Matrix
    Yint_fault = kron_reduce(Y_fault_phys)
    
    # 3. Post-Fault (Line Tripped)
    # Physical: Start with Loaded Pre-Fault, remove line
    Y_post1 = Y_pre_l.copy()
    Y_post1[idx_f, idx_f] -= y_orig_total
    Y_post1[idx_t, idx_t] -= y_orig_total
    Y_post1[idx_f, idx_t] = 0.0
    Y_post1[idx_t, idx_f] = 0.0
    
    # Reduced Post-Fault
    Yint_post = kron_reduce(Y_post1)
    
    # SAVE
    mdict = {
        "Yi_pre": Y_pre,
        "Y_fault": Y_fault_phys,   # 50% Split Topology
        "Y_post": Y_post1,         # Line Removed Topology
        "Yint_pre": Yint_pre,
        "Yint_fault": Yint_fault,
        "Yint_post": Yint_post
    }
    sio.savemat(mat_file, mdict)
    
    sio.savemat(mat_info, {
        "Fault_Type": "Line_50_Percent",
        "Fault_Location": f"Line_{trip_from}_{trip_to}",
        "Tripped_Line": [trip_from, trip_to],
        "Fault_Impedance": 0
    })
    print(f"Saved 50% Line Fault data to {mat_file}")